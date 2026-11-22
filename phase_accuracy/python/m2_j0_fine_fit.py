"""
m2_j0_fine_fit.py
==================
Locate the first J0 Bessel null (theoretical kp = 2.40483 rad) from the
11-point fine sweep, by fitting a parabola to the log-magnitude of the
carrier component.

Usage:
    python m2_j0_fine_fit.py --ref ref_carrier_kp0.csv

Expects files named m2_j0_fine_<kp>.csv matching the kp_list you flashed.
Edit KP_LIST below to match exactly what you put in firmware.
"""

import argparse
import numpy as np
from measure_fc_ref import load_tbs_csv   # reuse your existing loader

# EDIT to match exactly what you flashed (kp_min +/- 0.10, step 0.02)
KP_LIST = [2.30, 2.32, 2.34, 2.36, 2.38, 2.40, 2.42, 2.44, 2.46, 2.48, 2.50]

FC_NOMINAL = 5000.0
J0_NULL_THEORY = 2.40483


def carrier_amplitude(t, x, f_carrier=FC_NOMINAL):
    xa = x - np.mean(x)
    N = len(xa)
    dt = (t[-1] - t[0]) / (N - 1)
    w = np.blackman(N)
    X = np.abs(np.fft.rfft(xa * w))
    f = np.fft.rfftfreq(N, dt)
    k = int(np.argmin(np.abs(f - f_carrier)))
    return float(X[max(k - 2, 0):k + 3].max())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True,
                    help="ref_carrier_kp0.csv, for the 0 dB reference level")
    ap.add_argument("--pattern", default="m2_j0_fine_{kp:.2f}.csv")
    args = ap.parse_args()

    # unmodulated reference amplitude, for the suppression figure in dB
    t_ref, x_ref = load_tbs_csv(args.ref, channel=1)
    amp_ref = carrier_amplitude(t_ref, x_ref)

    kps, amps = [], []
    for kp in KP_LIST:
        fname = args.pattern.format(kp=kp)
        t, x = load_tbs_csv(fname, channel=1)
        a = carrier_amplitude(t, x)
        kps.append(kp)
        amps.append(a)
        print(f"  kp = {kp:.2f}   amp = {a:8.4f}   "
              f"{20*np.log10(a/amp_ref):+6.1f} dB")

    kps = np.asarray(kps)
    log_amp = np.log10(np.asarray(amps))

    # ---- parabolic fit across all 11 points ----
    # log_amp(kp) ~ a*kp^2 + b*kp + c  ->  vertex at kp = -b / (2a)
    a, b, c = np.polyfit(kps, log_amp, 2)
    if a <= 0:
        print("\nWARNING: fitted parabola opens downward -- the sweep did "
              "not bracket a minimum. Check the coarse-pass result and the "
              "kp_list actually flashed.")
    kp_null = -b / (2 * a)

    suppression_db = 20 * (np.polyval([a, b, c], kp_null) - np.log10(amp_ref))

    gain_error_pct = (kp_null - J0_NULL_THEORY) / J0_NULL_THEORY * 100.0

    print(f"\nFitted null location : kp = {kp_null:.4f} rad")
    print(f"Theoretical (J0)      : kp = {J0_NULL_THEORY:.5f} rad")
    print(f"Gain error             : {gain_error_pct:+.2f} %")
    print(f"Suppression at null    : {suppression_db:+.1f} dB "
          "(expect 30-40 dB)")

    if not (kps.min() < kp_null < kps.max()):
        print("\nWARNING: fitted minimum falls outside the swept range. "
              "The true null is probably beyond your fine sweep -- widen "
              "it and recapture rather than trusting this extrapolation.")

    if abs(suppression_db) < 20:
        print("\nWARNING: suppression under 20 dB is a shallow dip for a "
              "coarse-pass result. Re-check the coarse minimum before "
              "trusting this fit.")


if __name__ == "__main__":
    main()
