#!/usr/bin/env python3
"""
figure_7_11_b99_carson.py
---------------------------
Builds Figure 7.11: measured occupied bandwidth B99% compared against
Carson's rule, for the three phase-deviation test cases S3, S4 and S5.

    Carson's rule for phase modulation:  B = 2(kp + 1) * fm

Numbers are recomputed directly from the CSV data via occupied_bandwidth()
in analyse_spectrum.py, rather than hardcoded, so the figure always matches
whatever the main analysis currently reports.

The +-20% tolerance band used in the cross-check table (Section 7 of the
Python analysis procedure) is drawn as error bars on the Carson prediction,
so pass/fail is visible directly rather than requiring the reader to check
a percentage against a threshold stated only in text.

Requires analyse_spectrum.py and the same CSV files in the same folder.
"""

import numpy as np
import matplotlib.pyplot as plt

from analyse_spectrum import (
    collect_case, CASES, occupied_bandwidth, B99_LOW, B99_HIGH, FM,
)

CASES_TO_PLOT = ["S3", "S4", "S5"]
TOLERANCE = 0.20   # +-20%, matching the cross-check table


def carson_bandwidth(kp, fm=FM):
    return 2.0 * (kp + 1.0) * fm


def main():
    labels, measured, carson = [], [], []

    for case_id in CASES_TO_PLOT:
        cfg = CASES[case_id]
        res = collect_case(case_id, cfg)
        if res is None:
            print(f"{case_id}: no data found - skipping")
            continue

        f, Pxx = res["f"], res["Pxx"]
        b99, f_lo, f_hi = occupied_bandwidth(f, Pxx, B99_LOW, B99_HIGH)
        kp = cfg["kp"]
        c = carson_bandwidth(kp)
        dev = 100.0 * (b99 - c) / c

        labels.append(f"$k_p$ = {kp:.4f}")
        measured.append(b99)
        carson.append(c)

        flag = "PASS" if abs(dev) <= TOLERANCE * 100 else "FAIL"
        print(f"{case_id}: kp={kp:.4f}  measured={b99:.0f} Hz  "
              f"Carson={c:.0f} Hz  ({dev:+.1f}%)  [{flag}]")

    if not measured:
        print("Nothing to plot.")
        return

    x = np.arange(len(labels))
    width = 0.35
    carson_err = [c * TOLERANCE for c in carson]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(x - width / 2, measured, width,
           label="Measured $B_{99\\%}$", color="tab:blue")
    ax.bar(x + width / 2, carson, width, yerr=carson_err, capsize=6,
           label="Carson's rule ($\\pm$20% tolerance)", color="tab:orange")

    for xi, m, c in zip(x, measured, carson):
        dev = 100.0 * (m - c) / c
        y_top = max(m, c * (1 + TOLERANCE)) + 250
        ax.annotate(f"{dev:+.1f}%", (xi, y_top), ha="center", fontsize=9)

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Occupied bandwidth $B_{99\\%}$ (Hz)")
    ax.set_title("Figure 7.11 — Occupied bandwidth: measured vs. Carson's rule")
    ax.grid(True, alpha=0.3, axis="y")
    ax.legend()

    plt.tight_layout()
    plt.savefig("FIG_7_11_b99_carson.png", dpi=300)
    plt.savefig("FIG_7_11_b99_carson.pdf")
    print("\nSaved FIG_7_11_b99_carson.(png|pdf)")
    plt.show()


if __name__ == "__main__":
    main()
