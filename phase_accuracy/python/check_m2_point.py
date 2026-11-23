"""
check_m2_point.py
===================
Validates a single M2 (Bessel null) capture: checks capture mechanics
(timebase, clipping) and reports the carrier amplitude at fc via FFT --
the actual M2 measurand -- relative to the kp=0 reference.

diagnose.py's spectral checks are NOT reliable for M2 files, for two
reasons specific to this method: (1) at fm = 1000 Hz and kp >= 2.2, the
instantaneous frequency sweeps well outside diagnose.py's default
+/-20% band around fc, so legitimate PM signal shows reduced "power
near fc" even when nothing is wrong; (2) near an actual Bessel null the
carrier is INTENTIONALLY suppressed by 30-40 dB by design -- a low
carrier reading close to a null means the measurement is working, not
that the capture is bad.

Usage:
    python check_m2_point.py m2_j0_coarse_2.20.csv --kp 2.20 --ref ref_carrier_kp0.csv
"""

import argparse
import numpy as np
from measure_fc_ref import load_tbs_csv


def carrier_amplitude(t, x, f_carrier=5000.0):
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
    ap.add_argument("csv")
    ap.add_argument("--kp", type=float, required=True)
    ap.add_argument("--fc", type=float, default=5000.0)
    ap.add_argument("--ref", required=True,
                     help="ref_carrier_kp0.csv, for the dB reference level")
    args = ap.parse_args()

    t, x = load_tbs_csv(args.csv, channel=1)
    dt = (t[-1] - t[0]) / (len(t) - 1)

    print(f"record length  : {len(t)}  (expect 2000)")
    print(f"sample interval : {dt*1e6:.4f} us  (expect 8.0000 us)")

    # time-domain amplitude/clipping is still meaningful for M2 -- total
    # signal energy stays roughly constant as kp changes, since energy
    # moves into sidebands rather than disappearing
    xa = x - np.median(x)
    amp = (np.percentile(xa, 99.5) - np.percentile(xa, 0.5)) / 2.0
    top = np.percentile(xa, 99.9)
    bot = np.percentile(xa, 0.1)
    flat_top = np.mean(xa > 0.995 * top)
    flat_bot = np.mean(xa < 0.995 * bot)
    print(f"time-domain amplitude: {amp:.4f}")
    print(f"samples at extremes   : top {flat_top*100:.2f}%, "
          f"bottom {flat_bot*100:.2f}%")

    warn = []
    if flat_top > 0.06 or flat_bot > 0.06:
        warn.append("possible clipping -- check vertical scale")
    if amp < 0.3:
        warn.append("time-domain amplitude looks very low -- check DAC/"
                     "probe. Total signal energy should stay roughly "
                     "constant regardless of kp, unlike the carrier bin.")

    t_ref, x_ref = load_tbs_csv(args.ref, channel=1)
    amp_ref = carrier_amplitude(t_ref, x_ref, args.fc)
    amp_carrier = carrier_amplitude(t, x, args.fc)
    db = 20 * np.log10(amp_carrier / amp_ref)

    print(f"\ncarrier amplitude @ {args.fc:.0f} Hz : {amp_carrier:.4f}")
    print(f"suppression vs kp=0 reference     : {db:+.1f} dB")

    print("\n=== Verdict ===")
    if warn:
        for w in warn:
            print(f"  WARNING: {w}")
    else:
        print("PASS -- capture mechanics look fine.")

    print("\nNote: diagnose.py's 'power near fc' and slew-rate checks are")
    print("not meaningful for M2 (wide FM sweep at fm=1000Hz, and")
    print("intentional carrier suppression near the null) -- disregard")
    print("those specific warnings if you ran it on this file already.")


if __name__ == "__main__":
    main()
