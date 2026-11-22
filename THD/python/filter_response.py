#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
 FILTER RESPONSE  --  Figure 6, and the correction factors for case T4
===============================================================================

 Measures the real Sallen-Key low-pass response from a function generator
 sweep, fits the cutoff frequency and Q, and prints the de-embedding factors
 to use for the 5 kHz THD case.

 Input file: filter_sweep.csv, three columns, one header row

     freq_Hz,vin_pp,vout_pp
     100,1.004,1.002
     200,1.004,1.001
     ...

 Record Vpp in volts from the oscilloscope's automatic measurements, with
 CH1 on the filter input and CH2 on the filter output.

 Output:
     figures/fig6_filter_response.pdf and .png
     printed fitted fc and Q
     printed correction factors for harmonics 2 to 7 at f0 = 1 kHz and 5 kHz

===============================================================================
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# =============================================================================
# CONFIGURATION
# =============================================================================

SWEEP_FILE = "filter_sweep.csv"
FIG_DIR    = "figures"

FC_NOMINAL = 20100.0      # design value from R and C, Hz
Q_NOMINAL  = 0.7071       # Butterworth, 2nd order

F0_CASES   = [1000.0, 5000.0]   # the two THD test frequencies
N_HARM     = 7


# =============================================================================
# MODEL
# =============================================================================

def lowpass2(f, fc, q, g0):
    """
    Magnitude response of a general second-order low-pass.

        |H(f)| = g0 / sqrt( (1 - (f/fc)^2)^2 + (f / (q*fc))^2 )

    Both fc and q are fitted rather than assumed.  Real capacitors are
    typically 5 or 10 percent parts, so the built filter will not sit exactly
    on the Butterworth design values, and Q controls how much peaking or
    droop appears just below the corner.

    g0 is the passband gain.  It should come out near 1.000 for a unity-gain
    Sallen-Key stage; a fitted value far from 1 means a wiring or probe
    problem rather than a filter problem.
    """
    x = f / fc
    return g0 / np.sqrt((1.0 - x ** 2) ** 2 + (x / q) ** 2)


def db(x):
    return 20.0 * np.log10(np.asarray(x, dtype=float))


# =============================================================================
# MAIN
# =============================================================================

def main():
    if not os.path.exists(SWEEP_FILE):
        print(f"Cannot find {SWEEP_FILE}.")
        print("Create it with three columns: freq_Hz, vin_pp, vout_pp")
        return

    data = np.genfromtxt(SWEEP_FILE, delimiter=",", names=True)
    f = np.atleast_1d(data["freq_Hz"]).astype(float)
    vin = np.atleast_1d(data["vin_pp"]).astype(float)
    vout = np.atleast_1d(data["vout_pp"]).astype(float)

    order = np.argsort(f)
    f, vin, vout = f[order], vin[order], vout[order]

    gain = vout / vin
    gain_db = db(gain)

    print("=" * 70)
    print(" FILTER RESPONSE")
    print("=" * 70)
    print(f" {len(f)} points from {f[0]:.0f} Hz to {f[-1]/1000:.1f} kHz")

    # ---- fit ---------------------------------------------------------------
    # Weight by 1/gain so the stopband points, where the output is small,
    # still influence the fit.  Without weighting the fit is dominated by the
    # flat passband and the corner is poorly determined.
    sigma = 1.0 / np.maximum(gain, 1e-3)
    p0 = [FC_NOMINAL, Q_NOMINAL, 1.0]
    try:
        popt, pcov = curve_fit(lowpass2, f, gain, p0=p0, sigma=sigma,
                               maxfev=20000)
        perr = np.sqrt(np.diag(pcov))
    except Exception as exc:
        print(f" Fit failed: {exc}")
        return

    fc_fit, q_fit, g0_fit = popt
    fc_err, q_err, g0_err = perr

    print("\n Fitted parameters")
    print(f"   cutoff fc : {fc_fit/1000:.3f} kHz  +/- {fc_err/1000:.3f}")
    print(f"   Q         : {q_fit:.4f}  +/- {q_err:.4f}")
    print(f"   passband gain : {g0_fit:.4f}  +/- {g0_err:.4f}")

    print("\n Compared with the design values")
    print(f"   fc : nominal {FC_NOMINAL/1000:.2f} kHz, "
          f"measured {fc_fit/1000:.2f} kHz "
          f"({100*(fc_fit-FC_NOMINAL)/FC_NOMINAL:+.1f} %)")
    print(f"   Q  : nominal {Q_NOMINAL:.4f}, measured {q_fit:.4f} "
          f"({100*(q_fit-Q_NOMINAL)/Q_NOMINAL:+.1f} %)")

    if abs(g0_fit - 1.0) > 0.05:
        print("\n   WARNING: passband gain is more than 5 % from unity.")
        print("   Check the probe attenuation settings and the wiring before")
        print("   using this fit for anything.")

    resid = db(gain) - db(lowpass2(f, *popt))
    print(f"\n   fit residual : {np.std(resid):.3f} dB rms, "
          f"worst {np.max(np.abs(resid)):.3f} dB")

    # ---- correction factors -------------------------------------------------
    print("\n" + "=" * 70)
    print(" DE-EMBEDDING FACTORS FOR THE THD ANALYSIS")
    print("=" * 70)
    print(" Add these to each measured harmonic level in dBc before")
    print(" computing THD at Point B.  Positive numbers mean the filter")
    print(" removed distortion that the DAC really produced.\n")

    for f0 in F0_CASES:
        g1_fit = lowpass2(f0, *popt)
        g1_nom = lowpass2(f0, FC_NOMINAL, Q_NOMINAL, 1.0)
        print(f" f0 = {f0/1000:.0f} kHz")
        print(f"   {'harm':6}{'freq':>12}{'measured':>12}{'nominal':>12}"
              f"{'difference':>12}")
        for h in range(2, N_HARM + 1):
            fh = h * f0
            gm = lowpass2(fh, *popt)
            gn = lowpass2(fh, FC_NOMINAL, Q_NOMINAL, 1.0)
            cm = db(g1_fit / gm)
            cn = db(g1_nom / gn)
            print(f"   H{h:<5}{fh/1000:>10.0f} kHz{cm:>+11.2f} dB"
                  f"{cn:>+11.2f} dB{cm-cn:>+11.2f} dB")
        print()

    print(" If the measured and nominal columns differ by more than about")
    print(" 0.5 dB at any harmonic, use the measured filter in the THD")
    print(" analysis. Edit LPF_FCUT in thd_analysis.py to the fitted fc, and")
    print(" replace butter_gain with the fitted second-order model, since a")
    print(" Butterworth assumes Q = 0.707 exactly.\n")

    # ---- figure -------------------------------------------------------------
    os.makedirs(FIG_DIR, exist_ok=True)
    fsmooth = np.logspace(np.log10(f[0] * 0.8), np.log10(f[-1] * 1.2), 800)

    fig, ax = plt.subplots(figsize=(9, 5.5))

    ax.semilogx(fsmooth / 1000,
                db(lowpass2(fsmooth, FC_NOMINAL, Q_NOMINAL, 1.0)),
                ls="--", lw=1.2, color="0.55",
                label=f"design, $f_c$ = {FC_NOMINAL/1000:.1f} kHz, "
                      f"Q = {Q_NOMINAL:.3f}")

    ax.semilogx(fsmooth / 1000, db(lowpass2(fsmooth, *popt)),
                lw=1.6, color="#1f77b4",
                label=f"fitted, $f_c$ = {fc_fit/1000:.1f} kHz, "
                      f"Q = {q_fit:.3f}")

    ax.semilogx(f / 1000, gain_db, "o", ms=5, color="#d62728",
                label="measured", zorder=5)

    # mark where the harmonics land for each test frequency
    styles = {1000.0: ("#2ca02c", "-."), 5000.0: ("#ff7f0e", ":")}
    for f0 in F0_CASES:
        col, ls = styles.get(f0, ("0.3", ":"))
        for h in range(2, N_HARM + 1):
            fh = h * f0
            if fh < f[0] * 0.8 or fh > f[-1] * 1.2:
                continue
            ax.axvline(fh / 1000, color=col, ls=ls, lw=0.8, alpha=0.55)
        ax.plot([], [], color=col, ls=ls, lw=0.8,
                label=f"harmonics 2-7 of $f_0$ = {f0/1000:.0f} kHz")

    ax.axhline(-3.0, color="0.7", lw=0.8, ls=":")
    ax.text(f[0] / 1000, -3.0 + 0.4, " -3 dB", fontsize=8, color="0.5",
            va="bottom")

    ax.set_xlabel("Frequency (kHz)")
    ax.set_ylabel("Gain (dB)")
    ax.set_title("Measured response of the Sallen-Key second-order low-pass "
                 "filter,\nwith the harmonic positions of both THD test "
                 "frequencies")
    ax.set_ylim(min(-35, float(np.min(gain_db)) - 3), 5)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, loc="lower left")

    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG_DIR, f"fig6_filter_response.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f" Figure written to {FIG_DIR}/fig6_filter_response.pdf and .png")


if __name__ == "__main__":
    main()
