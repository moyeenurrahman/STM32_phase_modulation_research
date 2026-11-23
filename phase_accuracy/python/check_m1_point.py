"""
check_m1_point.py
===================
Analyzes a single M1 sweep capture: extracts the phase trajectory via
Hilbert baseband conversion, fits the group delay and constant offset,
and reports the measured beta, gain error, and RMS/peak phase error.

This is the actual per-point M1 result -- reuse it for every kp in the
sweep (0.25 through 4.00), all three repeats.

Usage:
    python check_m1_point.py m1_kp0p25_r1.csv --kp 0.25 --fc 5000.12
"""

import argparse
import numpy as np
from scipy.signal import hilbert, butter, filtfilt
from scipy.optimize import minimize_scalar
from measure_fc_ref import load_tbs_csv

FM_HZ = 500.0
BW_HZ = 3000.0            # post-detection bandwidth -- state with every result
TRIM_CYCLES = 10
NOISE_FLOOR_DEG = 0.30    # measured in R3, Average x16


def wrap(a):
    return (a + np.pi) % (2 * np.pi) - np.pi


def circular_mean(angles):
    return float(np.angle(np.mean(np.exp(1j * angles))))


def extract_phase(t, x, fc, bw_hz, trim_cycles=10):
    fs = 1.0 / ((t[-1] - t[0]) / (len(t) - 1))
    x_ac = x - np.mean(x)
    z = hilbert(x_ac)
    zb = z * np.exp(-1j * 2 * np.pi * fc * t)
    b, a = butter(4, bw_hz / (fs / 2.0))
    zb = filtfilt(b, a, zb.real) + 1j * filtfilt(b, a, zb.imag)
    theta = np.unwrap(np.angle(zb))
    env = np.abs(zb)
    A = np.vstack([t, np.ones_like(t)]).T
    coef, *_ = np.linalg.lstsq(A, theta, rcond=None)
    theta = theta - A @ coef
    n = int(trim_cycles * fs / fc)
    sl = slice(n, -n)
    return t[sl], theta[sl], env[sl]


def fit_delay_offset(t, theta_meas, ideal_fn, tau_bound=5e-4):
    def cost(tau):
        e = wrap(theta_meas - ideal_fn(t - tau))
        return float(np.mean(wrap(e - circular_mean(e)) ** 2))
    res = minimize_scalar(cost, bounds=(-tau_bound, tau_bound), method="bounded")
    tau = float(res.x)
    e = wrap(theta_meas - ideal_fn(t - tau))
    offset = circular_mean(e)
    err = wrap(e - offset)
    return tau, offset, err


def measure_beta(t, theta_meas, fm):
    s = np.sin(2 * np.pi * fm * t)
    c = np.cos(2 * np.pi * fm * t)
    A = np.vstack([s, c, np.ones_like(t)]).T
    coef, *_ = np.linalg.lstsq(A, theta_meas, rcond=None)
    return float(np.hypot(coef[0], coef[1]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--kp", type=float, required=True, help="commanded kp, rad")
    ap.add_argument("--fc", type=float, required=True, help="measured fc from R1")
    ap.add_argument("--bw", type=float, default=BW_HZ)
    args = ap.parse_args()

    t, x = load_tbs_csv(args.csv, channel=1)
    fs = 1.0 / ((t[-1] - t[0]) / (len(t) - 1))

    # envelope ripple, computed the same way as envelope_check.py
    xa = x - np.mean(x)
    zb_raw = hilbert(xa) * np.exp(-1j * 2 * np.pi * args.fc * t)
    b, a = butter(4, args.bw / (fs / 2.0))
    zb_f = filtfilt(b, a, zb_raw.real) + 1j * filtfilt(b, a, zb_raw.imag)
    n = int(TRIM_CYCLES * fs / args.fc)
    env = np.abs(zb_f)[n:-n]
    ripple_pct = float((env.max() - env.min()) / env.mean() * 100)

    t_p, theta, _ = extract_phase(t, x, args.fc, args.bw, TRIM_CYCLES)

    beta_hat = measure_beta(t_p, theta, FM_HZ)
    gain_error_pct = (beta_hat / args.kp - 1) * 100 if args.kp != 0 else float("nan")

    ideal = lambda tt: args.kp * np.sin(2 * np.pi * FM_HZ * tt)
    tau, offset, err = fit_delay_offset(t_p, theta, ideal)

    rms_deg = float(np.degrees(np.sqrt(np.mean(err ** 2))))
    peak_deg = float(np.degrees(np.max(np.abs(err))))

    print(f"=== {args.csv}  (kp = {args.kp:.2f} rad) ===")
    print(f"Envelope ripple    : {ripple_pct:.2f} %   (limit 5%)")
    print(f"Measured beta       : {beta_hat:.4f} rad")
    print(f"Gain error           : {gain_error_pct:+.2f} %")
    print(f"Fitted delay         : {tau*1e6:+.2f} us   (nuisance parameter, see note)")
    print(f"Fitted offset        : {np.degrees(offset):+.2f} deg   "
          f"(nuisance parameter, see note)")
    print(f"RMS phase error       : {rms_deg:.3f} deg   "
          f"(noise floor: {NOISE_FLOOR_DEG:.2f} deg)")
    print(f"Peak phase error      : {peak_deg:.3f} deg")

    warn = []
    if ripple_pct > 5:
        warn.append("envelope ripple exceeds 5% -- check triggering/averaging")

    print("\n=== Verdict ===")
    if warn:
        for w in warn:
            print(f"  WARNING: {w}")
    else:
        print("PASS -- looks like a good M1 data point.")

    print("\nNote: fitted delay/offset only align the trajectories before")
    print("computing RMS/peak error -- they are relative to this capture's")
    print("own trigger point, not an absolute time origin (M1 is self-")
    print("triggered, not synced to the DMA buffer wrap). They will not")
    print("match R2's filter-only delay/offset and should not be compared")
    print("to it. R2 remains the authoritative filter characterization.")

    if rms_deg < NOISE_FLOOR_DEG:
        print(f"\nNOTE: RMS error is at or below the measured noise floor -- "
              f"report this point as instrument-limited, quoting "
              f"{NOISE_FLOOR_DEG:.2f} deg, rather than the raw RMS value.")


if __name__ == "__main__":
    main()
