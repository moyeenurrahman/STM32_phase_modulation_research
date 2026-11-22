#!/usr/bin/env python3
"""
figure_7_10_bessel_fit.py
---------------------------
Builds Figure 7.10: measured PM sideband amplitudes at orders n = 0..6
compared against the Bessel function pattern |J_n(kp)|, for the three
phase-deviation test cases S3, S4 and S5.

The theoretical curve uses the RECOVERED kp (the least-squares fit from
recover_kp() in analyse_spectrum.py), not the commanded kp. This is what
makes the plot show the same thing the printed "fit residual" number
describes: how closely the measured amplitude PATTERN matches a Bessel
shape, independent of any small deviation-index error.

Both measured and theoretical arrays are normalised by their own L2 norm,
identically to how recover_kp() performs the fit internally, so the
figure and the printed residual are describing exactly the same quantity.

Plotted in dB rather than linear amplitude, because a linear scale would
show only the first two or three orders - everything beyond that is
visually indistinguishable from zero even where the fit is excellent to
within a fraction of a dB.

Requires analyse_spectrum.py and the same CSV files in the same folder.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv

from analyse_spectrum import (
    collect_case, CASES, lobe_amplitude, recover_kp,
    FC, FM, N_SIDEBANDS,
)

CASES_TO_PLOT = ["S3", "S4", "S5"]
DB_FLOOR = -60.0     # y-axis floor; anything below this is clipped for display


def get_amplitudes(f, Pxx):
    """
    Same extraction method used in analyse_spectrum.py's analyse_family_n:
    carrier amplitude at n=0, then the mean of the lower/upper sideband
    amplitude at each order n=1..N_SIDEBANDS.
    """
    amps = [lobe_amplitude(f, Pxx, FC)]
    for n in range(1, N_SIDEBANDS + 1):
        a_lo = lobe_amplitude(f, Pxx, FC - n * FM)
        a_hi = lobe_amplitude(f, Pxx, FC + n * FM)
        amps.append(0.5 * (a_lo + a_hi))
    return np.array(amps)


def to_db_clipped(x, floor_db=DB_FLOOR):
    x = np.maximum(x, 10 ** (floor_db / 20.0))
    return 20.0 * np.log10(x)


def main():
    orders = np.arange(0, N_SIDEBANDS + 1)
    fig, axes = plt.subplots(3, 1, figsize=(8, 11), sharex=True)

    for ax, case_id in zip(axes, CASES_TO_PLOT):
        cfg = CASES[case_id]
        res = collect_case(case_id, cfg)
        if res is None:
            ax.set_title(f"{case_id}: no data found")
            continue

        f, Pxx = res["f"], res["Pxx"]
        amps = get_amplitudes(f, Pxx)

        kp_fit, resid = recover_kp(amps)
        kp_cmd = cfg["kp"]

        a_norm = amps / np.linalg.norm(amps)
        theory = np.abs(jv(orders, kp_fit))
        theory_norm = theory / np.linalg.norm(theory)

        a_db = to_db_clipped(a_norm)
        theory_db = to_db_clipped(theory_norm)

        width = 0.35
        ax.bar(orders - width / 2, a_db, width,
               label="Measured", color="tab:blue")
        ax.bar(orders + width / 2, theory_db, width,
               label=f"Bessel theory (fitted $k_p$={kp_fit:.4f})",
               color="tab:orange", alpha=0.85)

        err_pct = 100.0 * (kp_fit - kp_cmd) / kp_cmd if kp_cmd else 0.0
        ax.set_title(
            f"{case_id}:  commanded $k_p$ = {kp_cmd:.4f} rad,  "
            f"recovered $k_p$ = {kp_fit:.4f} rad  "
            f"({err_pct:+.2f}%),  fit residual = {resid:.4f}"
        )
        ax.set_ylabel("Normalised amplitude (dB)")
        ax.set_ylim(DB_FLOOR, 5)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, loc="upper right")

        print(f"{case_id}: commanded={kp_cmd:.4f}  recovered={kp_fit:.4f}  "
              f"residual={resid:.4f}")

    axes[-1].set_xlabel("Sideband order  n")
    axes[-1].set_xticks(orders)

    fig.suptitle("Figure 7.10 — Measured sideband amplitudes vs. Bessel prediction",
                 y=1.0)
    plt.tight_layout()
    plt.savefig("FIG_7_10_bessel_fit.png", dpi=300, bbox_inches="tight")
    plt.savefig("FIG_7_10_bessel_fit.pdf", bbox_inches="tight")
    print("\nSaved FIG_7_10_bessel_fit.(png|pdf)")
    plt.show()


if __name__ == "__main__":
    main()
