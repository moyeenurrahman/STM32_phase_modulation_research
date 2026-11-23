"""
check_r2.py
============
Validates the R2 reference capture: measures the LPF's phase offset
relative to the raw DAC output (PA4), via a least-squares sine fit on
each channel at the measured carrier frequency.

At kp = 0 (unmodulated), a constant phase shift and an equivalent time
delay cannot be told apart from a single-tone measurement -- they are
the same quantity, expressed two ways. The independent fit that later
separates them (fit_delay_offset) needs modulated (kp != 0) data from
the M1 sweep. This script is a sanity check only: it confirms the
filter + DAC zero-order-hold behave roughly as expected before you
trust that later fit.

Usage:
    python check_r2.py TEK00000.CSV TEK00001.CSV --fc 5000.12

The first file must be CH1 (LPF output), the second CH2 (raw PA4).
"""

import argparse
import numpy as np
from measure_fc_ref import load_tbs_csv


def sine_phase(t, x, f):
    x = x - np.mean(x)
    A = np.vstack([np.sin(2 * np.pi * f * t), np.cos(2 * np.pi * f * t),
                   np.ones_like(t)]).T
    c, *_ = np.linalg.lstsq(A, x, rcond=None)
    amp = float(np.hypot(c[0], c[1]))
    phase = float(np.arctan2(c[1], c[0]))
    return amp, phase


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ch1_csv", help="LPF output (filtered)")
    ap.add_argument("ch2_csv", help="raw PA4 (unfiltered)")
    ap.add_argument("--fc", type=float, required=True)
    args = ap.parse_args()

    t1, x1 = load_tbs_csv(args.ch1_csv, channel=1)
    t2, x2 = load_tbs_csv(args.ch2_csv, channel=1)

    amp1, ph1 = sine_phase(t1, x1, args.fc)
    amp2, ph2 = sine_phase(t2, x2, args.fc)

    offset_deg = float(np.degrees(np.angle(np.exp(1j * (ph1 - ph2)))))
    equiv_delay_us = (offset_deg / 360.0) / args.fc * 1e6

    print(f"CH1 (LPF output) amplitude : {amp1:.4f}")
    print(f"CH2 (raw PA4)    amplitude : {amp2:.4f}")
    print(f"Constant phase offset       : {offset_deg:+.2f} deg")
    print(f"Equivalent time delay       : {equiv_delay_us:+.2f} us")
    print()
    print("Expected: the filter alone contributes about -20.7 deg; the DAC")
    print("zero-order-hold adds roughly another -9 deg, so a total near")
    print("-30 deg (range about -15 to -40 deg) is normal. This is a")
    print("sanity check only -- the definitive offset/delay split comes")
    print("later from the modulated (kp != 0) fit during the M1 sweep.")

    if not (-45 <= offset_deg <= -10):
        print("\nWARNING: offset is outside the expected range. Check "
              "filter component values (R1, R2, C1, C2) before trusting "
              "this reference.")
    else:
        print("\nPASS -- offset is in the expected range.")


if __name__ == "__main__":
    main()
