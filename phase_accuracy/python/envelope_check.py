"""
envelope_check.py
==================
R4 envelope diagnostic: computes the instantaneous envelope |zb(t)| via
Hilbert transform and complex-baseband conversion, and reports the ripple
percentage. Per the methodology, ripple over ~5% invalidates the phase
extraction and must be fixed before proceeding to the M1 sweep.

Usage:
    python envelope_check.py TEK00000.CSV --fc 5000.12
"""

import argparse
import numpy as np
from scipy.signal import hilbert, butter, filtfilt
from measure_fc_ref import load_tbs_csv

FM_HZ = 500.0
TRIM_CYCLES = 10


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--fc", type=float, required=True,
                    help="measured carrier freq from R1 (fc_ref.json)")
    ap.add_argument("--bw", type=float, default=3000.0,
                    help="post-detection bandwidth, Hz")
    args = ap.parse_args()

    t, x = load_tbs_csv(args.csv, channel=1)
    fs = 1.0 / ((t[-1] - t[0]) / (len(t) - 1))

    xa = x - np.mean(x)
    z = hilbert(xa)
    zb = z * np.exp(-1j * 2 * np.pi * args.fc * t)

    b, a = butter(4, args.bw / (fs / 2.0))
    zb = filtfilt(b, a, zb.real) + 1j * filtfilt(b, a, zb.imag)

    n = int(TRIM_CYCLES * fs / args.fc)
    env = np.abs(zb)[n:-n]

    ripple_pct = (env.max() - env.min()) / env.mean() * 100.0

    print(f"Envelope mean   : {env.mean():.4f}")
    print(f"Envelope min/max: {env.min():.4f} / {env.max():.4f}")
    print(f"Ripple          : {ripple_pct:.2f} %  (limit: 5 %)")

    if ripple_pct > 5.0:
        print("\nFAIL -- envelope is not flat enough. Do not trust the phase "
              "extraction from this capture. Improve triggering before "
              "recapturing (wider sync pulse, check trigger level/holdoff).")
    else:
        print("\nPASS -- envelope is flat enough for phase extraction.")


if __name__ == "__main__":
    main()
