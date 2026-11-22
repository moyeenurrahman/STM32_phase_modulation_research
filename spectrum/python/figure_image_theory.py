#!/usr/bin/env python3
"""
figure_7_6_image_theory.py
---------------------------
Builds Figure 7.6: measured image attenuation at Test Point B (case S2)
compared against two theoretical predictions.

    1. Ideal:        sinc rolloff x 2nd-order Butterworth response.
                      Valid if the op-amp has enough loop gain to realise
                      the filter exactly, at every frequency plotted.

    2. Feedthrough:   at 95/105 kHz the MCP6002 (1 MHz GBW) has almost no
                      loop gain left. C1 then couples the input side
                      through to the output directly, and the stopband
                      floors at the R1/Xc1 divider ratio instead of
                      continuing to fall with the ideal Butterworth curve.

The feedthrough model is plotted ONLY at 95 kHz and 105 kHz, as discrete
markers, not as a continuous curve. It was derived and checked against
measurement at exactly those two frequencies (S1 vs S2); extending the
same formula down toward the passband would not be justified, since nothing
here confirms the op-amp has lost gain at, say, 30 kHz.

Requires analyse_spectrum.py and the same CSV files in the same folder.
"""

import numpy as np
import matplotlib.pyplot as plt

from analyse_spectrum import collect_case, CASES, to_db

# ---- filter component values (from the methodology chapter) ----
R1 = 5600.0            # ohm
C1 = 2e-9               # farad
FC_FILTER = 20100.0     # Hz, measured -3 dB point
FS_DAC = 100000.0       # Hz, DAC sample rate

F_PLOT_LOW = 20_000.0
F_PLOT_HIGH = 115_000.0


def sinc_db(f):
    """ZOH sinc roll-off relative to DC, in dB."""
    x = f / FS_DAC
    return 20.0 * np.log10(np.abs(np.sinc(x)) + 1e-12)


def butterworth_ideal_db(f, fc=FC_FILTER, order=2):
    """Ideal maximally-flat-magnitude response, in dB."""
    return -10.0 * np.log10(1.0 + (f / fc) ** (2 * order))


def feedthrough_db(f, R1=R1, C1=C1):
    """
    C1 feedthrough estimate. Once the op-amp has no usable loop gain,
    C1 couples node A almost directly to the output, forming a simple
    R1-C1 divider. Only meaningful well beyond the op-amp's usable
    bandwidth - do not extrapolate this toward the passband.
    """
    Xc1 = 1.0 / (2.0 * np.pi * f * C1)
    return 20.0 * np.log10(Xc1 / R1)


def main():
    res = collect_case("S2", CASES["S2"])
    if res is None:
        print("Could not load S2 - check DATA_DIR and CASES in analyse_spectrum.py")
        return

    f, Pxx = res["f"], res["Pxx"]
    Pxx_dB = to_db(Pxx)

    # Normalise the measured trace to its own peak (the carrier bin).
    # This is a visualisation choice: it will not reproduce the exact
    # dBc figures in the results table, which were computed from summed
    # lobe power rather than a single bin, but the two agree to within a
    # few tenths of a dB. What matters for this figure is the SHAPE of
    # the curve near 90-110 kHz relative to the two theoretical lines,
    # not an exact absolute number - the exact numbers are already in
    # the results text and table.
    measured_dbc = Pxx_dB - np.max(Pxx_dB)

    m = (f >= F_PLOT_LOW) & (f <= F_PLOT_HIGH)
    f_plot = f[m]

    ideal_curve = sinc_db(f_plot) + butterworth_ideal_db(f_plot)

    f_images = np.array([95_000.0, 105_000.0])
    feedthrough_curve = sinc_db(f_images) + feedthrough_db(f_images)

    print(f"Used {res['n_used']} averaged records for S2 "
          f"(found {res['n_files']}, {res['n_dup']} duplicates)")
    print(f"\nIdeal prediction   at 95/105 kHz: "
          f"{ideal_curve[np.argmin(np.abs(f_plot-95000))]:.1f} / "
          f"{ideal_curve[np.argmin(np.abs(f_plot-105000))]:.1f} dB")
    print(f"Feedthrough model  at 95/105 kHz: "
          f"{feedthrough_curve[0]:.1f} / {feedthrough_curve[1]:.1f} dB")

    plt.figure(figsize=(9, 5))
    plt.plot(f_plot / 1e3, measured_dbc[m], linewidth=1.0,
             label="Measured (Point B, S2)", color="tab:blue")
    plt.plot(f_plot / 1e3, ideal_curve, linewidth=1.2, linestyle="--",
             label="Ideal: sinc x 2nd-order Butterworth", color="tab:green")
    plt.scatter(f_images / 1e3, feedthrough_curve, marker="x", s=80,
                color="tab:red", zorder=5,
                label="C1 feedthrough estimate (95/105 kHz only)")

    for fx, fy in zip(f_images, feedthrough_curve):
        plt.annotate(f"{fy:.1f} dB", (fx / 1e3, fy),
                     textcoords="offset points", xytext=(6, -10),
                     fontsize=8, color="tab:red")

    plt.axvline(95, color="grey", linestyle=":", alpha=0.5)
    plt.axvline(105, color="grey", linestyle=":", alpha=0.5)
    plt.xlabel("Frequency (kHz)")
    plt.ylabel("Level relative to carrier (dB)")
    plt.title("Figure 7.6 — Measured vs. predicted image attenuation")
    plt.xlim(F_PLOT_LOW / 1e3, F_PLOT_HIGH / 1e3)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig("FIG_7_6_image_theory.png", dpi=300)
    plt.savefig("FIG_7_6_image_theory.pdf")
    print("\nSaved FIG_7_6_image_theory.(png|pdf)")
    plt.show()


if __name__ == "__main__":
    main()
