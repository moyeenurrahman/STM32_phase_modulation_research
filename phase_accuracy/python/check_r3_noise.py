"""
check_r3_noise.py
===================
Validates the R3 noise-floor capture (DAC held static at MID_CODE, no
DMA, no timer) and computes the resulting phase measurement floor, per
Section 9 of the methodology: sigma_phi ~= sigma_n / A.

This capture should look like flat noise with no real 5 kHz content --
diagnose.py's spectral checks are not meaningful here and will likely
(correctly) flag it as "not a sine". Ignore that. This script does the
right checks for a noise-floor file instead.

Usage:
    python check_r3_noise.py R3_CH1.CSV --codes-per-volt 1237.4

--codes-per-volt: DAC codes per volt, used to convert the measured
    noise from volts into codes so it can be compared against the
    methodology's code-based floor derivation. Calibrate this from your
    own session data: AMP_CODES / measured_amplitude_volts (e.g.
    1200 / 0.9698 = 1237.4 from your R2 result). Falls back to the
    nominal 4095/3.3 = 1240.9 if not given.
"""

import argparse
import numpy as np
from measure_fc_ref import load_tbs_csv

AMP_CODES = 1200          # peak amplitude used for all phase captures
NOMINAL_CODES_PER_VOLT = 4095 / 3.3


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--codes-per-volt", type=float,
                    default=NOMINAL_CODES_PER_VOLT)
    args = ap.parse_args()

    t, x = load_tbs_csv(args.csv, channel=1)
    dt = (t[-1] - t[0]) / (len(t) - 1)

    print(f"record length      : {len(t)}  (expect 2000)")
    print(f"sample interval     : {dt*1e6:.4f} us  (expect 8.0000 us)")

    # this is genuinely flat noise -- RMS of the fluctuation, not a sine fit
    x_ac = x - np.mean(x)
    rms_v = float(np.std(x_ac))
    pk_v = float(np.max(np.abs(x_ac)))

    rms_codes = rms_v * args.codes_per_volt
    pk_codes = pk_v * args.codes_per_volt

    print(f"\nnoise RMS           : {rms_v*1e3:.3f} mV  ({rms_codes:.3f} codes)")
    print(f"noise peak          : {pk_v*1e3:.3f} mV  ({pk_codes:.3f} codes)")

    # check for a dominant tone rather than true broadband noise -- mains
    # hum, digital switching pickup, etc. A real noise floor should NOT
    # have one frequency carrying most of the power.
    N = len(x_ac)
    w = np.blackman(N)
    X = np.abs(np.fft.rfft(x_ac * w))
    f = np.fft.rfftfreq(N, dt)
    k_peak = int(np.argmax(X[1:]) + 1)     # skip DC bin
    frac_at_peak = float(X[k_peak] ** 2 / np.sum(X[1:] ** 2))

    print(f"\nstrongest tone      : {f[k_peak]:.0f} Hz, "
          f"{frac_at_peak*100:.1f}% of AC power")
    if frac_at_peak > 0.15:
        print("WARNING: one frequency dominates the noise -- this looks "
              "like pickup or interference, not a clean noise floor. "
              "Check grounding and nearby switching sources before "
              "trusting this floor.")

    # Section 9 phase floor: sigma_phi ~= sigma_n / A, A in codes
    sigma_phi_rad = rms_codes / AMP_CODES
    sigma_phi_deg = np.degrees(sigma_phi_rad)

    print(f"\nPhase floor (this capture, no averaging):")
    print(f"  sigma_phi = {sigma_phi_deg:.4f} deg")
    print(f"  (methodology's standard-8-bit reference figure: ~0.16 deg)")

    print(f"\nExpected floor after averaging (reduces by sqrt(N)):")
    for n in (16, 64):
        print(f"  Average x{n:<3}: {sigma_phi_deg/np.sqrt(n):.4f} deg")

    print("\nNote: diagnose.py's spectral/amplitude checks are not "
          "meaningful for this file and may show warnings -- that is "
          "expected for a noise-floor capture, not a fault.")


if __name__ == "__main__":
    main()
