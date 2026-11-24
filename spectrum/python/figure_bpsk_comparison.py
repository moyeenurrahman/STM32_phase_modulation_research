#!/usr/bin/env python3
"""
figure_7_14_bpsk_comparison.py
--------------------------------
Builds Figure 7.14: BPSK spectrum comparison between abrupt phase
switching (S6N) and the ramped transition (S7N).

This is a two-panel figure, not a single overlay, because the two cases
tell a mixed story that a single zoomed plot would hide:

    Panel A (main lobe, zoomed): shows the null-to-null width. S7N's
    nulls are measurably WIDER than S6N's, not narrower - the ramp
    introduces a transient frequency offset during each phase transition,
    which broadens the main lobe rather than containing it.

    Panel B (wide view, harmonics marked): shows H2 through H6. S7N's H2
    is worse than S6N's, not better.

Report this as the negative/mixed result it is. Do not caption this
figure as "the ramp reduces spectral spread" - the data does not support
that conclusion for phase-domain ramping at these settings.

Requires analyse_spectrum.py and the same CSV files in the same folder.
"""

import os

import numpy as np
import matplotlib.pyplot as plt

from analyse_spectrum import (
    OUT_DIR,
    collect_case, CASES, to_db, find_null, lobe_power,
    FC, RS, N_HARMONICS,
)

CASES_TO_PLOT = ["S6N", "S7N"]
COLORS = {"S6N": "tab:blue", "S7N": "tab:orange"}
LABELS = {"S6N": "Abrupt switching (S6N)", "S7N": "Phase ramp (S7N)"}


def main():
    results = {}
    for case_id in CASES_TO_PLOT:
        res = collect_case(case_id, CASES[case_id])
        if res is None:
            print(f"{case_id}: no data found")
            continue
        results[case_id] = res

    if len(results) < 2:
        print("Need both S6N and S7N to build the comparison.")
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    # ---- Panel A: main lobe, zoomed, nulls marked ----
    print("Main lobe:")
    for case_id, res in results.items():
        f, P = res["f"], res["Pxx"]
        Pdb = to_db(P)
        ax1.plot(f / 1e3, Pdb, linewidth=1.1, color=COLORS[case_id],
                 label=LABELS[case_id])

        f_nl, _ = find_null(f, P, FC - RS)
        f_nu, _ = find_null(f, P, FC + RS)
        ax1.axvline(f_nl / 1e3, color=COLORS[case_id], linestyle=":", alpha=0.5)
        ax1.axvline(f_nu / 1e3, color=COLORS[case_id], linestyle=":", alpha=0.5)

        print(f"  {case_id}: null-to-null = {f_nu - f_nl:.0f} Hz "
              f"({f_nl:.0f} - {f_nu:.0f} Hz)")

    ax1.set_xlim(0, 10)
    ax1.set_xlabel("Frequency (kHz)")
    ax1.set_ylabel("PSD (dB(V$^2$/Hz))")
    ax1.set_title("Main lobe — dotted lines mark measured nulls")
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=8)

    # ---- Panel B: wide view, harmonics marked ----
    print("\nHarmonics (relative to carrier bin):")
    for case_id, res in results.items():
        f, P = res["f"], res["Pxx"]
        Pdb = to_db(P)
        ax2.plot(f / 1e3, Pdb, linewidth=1.1, color=COLORS[case_id],
                 label=LABELS[case_id])

        p_carrier = lobe_power(f, P, FC)
        p_h2 = lobe_power(f, P, 2 * FC)
        print(f"  {case_id}: H2 = {to_db(p_h2):.2f} dB(V^2)  "
              f"({to_db(p_h2) - to_db(p_carrier):+.2f} dBc)")

    for n in range(2, N_HARMONICS + 1):
        ax2.axvline(n * FC / 1e3, color="grey", linestyle=":", alpha=0.4)
        ax2.text(n * FC / 1e3, 0.97, f"H{n}",
                  transform=ax2.get_xaxis_transform(),
                  ha="center", va="top", fontsize=7, color="grey")

    ax2.set_xlim(0, 32)
    ax2.set_xlabel("Frequency (kHz)")
    ax2.set_ylabel("PSD (dB(V$^2$/Hz))")
    ax2.set_title("Wide view — harmonics marked")
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=8)

    fig.suptitle(
        "Figure 7.14 — BPSK spectrum: abrupt vs. ramped phase transition",
        y=1.02,
    )
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "bpsk_comparison.png"), dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(OUT_DIR, "bpsk_comparison.pdf"), bbox_inches="tight")
    print("\nSaved bpsk_comparison.(png|pdf) in reproduced_plots/")
    plt.show()


if __name__ == "__main__":
    main()
