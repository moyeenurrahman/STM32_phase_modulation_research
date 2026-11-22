"""
quick_coarse_check.py
======================
Quick look at the M2 J0 coarse sweep: prints carrier amplitude and
suppression (dB relative to the kp=0 reference) for each coarse-pass
capture, so you can read off which kp gave the deepest dip before
running the fine sweep.

Requires measure_fc_ref.py in the same folder (for the CSV loader).
Run with the working directory set to the folder containing your
m2_j0_coarse_*.csv files and ref_carrier_kp0.csv.
"""

import numpy as np
from measure_fc_ref import load_tbs_csv

KP_VALS = [2.20, 2.25, 2.30, 2.35, 2.40, 2.45, 2.50, 2.55, 2.60]
FILE_PATTERN = "m2_j0_coarse_{kp:.2f}.csv"
REF_FILE = "ref_carrier_kp0.csv"          # kp = 0 reference, from R1


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
    t_ref, x_ref = load_tbs_csv(REF_FILE, channel=1)
    amp_ref = carrier_amplitude(t_ref, x_ref)
    print(f"Reference (kp=0) amplitude: {amp_ref:.4f}\n")

    results = []
    for kp in KP_VALS:
        fname = FILE_PATTERN.format(kp=kp)
        t, x = load_tbs_csv(fname, channel=1)
        amp = carrier_amplitude(t, x)
        db = 20 * np.log10(amp / amp_ref)
        results.append((kp, amp, db))
        print(f"kp = {kp:.2f}   amp = {amp:8.4f}   {db:+6.1f} dB")

    kp_min, amp_min, db_min = min(results, key=lambda r: r[1])
    print(f"\nCoarse minimum: kp = {kp_min:.2f}  "
          f"(amp = {amp_min:.4f}, {db_min:+.1f} dB)")
    print("Center the fine sweep (Section 7.2 step 5) on this value.")


if __name__ == "__main__":
    main()
