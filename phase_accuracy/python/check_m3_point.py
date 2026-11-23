"""
check_m3_point.py
===================
Analyzes an M3 (BPSK) capture: locates t_sync from the CH2 sync pulse
(the pulse's rising edge is the true DMA buffer-wrap moment, found
directly from the data -- independent of wherever the scope's own
trigger reference happens to sit), extracts the phase trajectory from
CH1, compares it against the ideal bit pattern, and reports state
separation and steady-state RMS phase error.

Usage:
    python check_m3_point.py m3_ch1.csv m3_ch2.csv --fc 5000.12

CH1 = LPF output, CH2 = PB0 sync pulse. M3 is the one method where
CH2's actual data is required, not just its role as a trigger.
"""

import argparse
import numpy as np
from scipy.signal import hilbert, butter, filtfilt
from measure_fc_ref import load_tbs_csv

RS = 250.0
GUARD = 400e-6
BW_HZ = 2500.0
TRIM_CYCLES = 10
bits = [0, 1, 0, 1, 1, 0, 1, 0, 0, 1]


def wrap(a):
    return (a + np.pi) % (2 * np.pi) - np.pi


def circular_mean(angles):
    return float(np.angle(np.mean(np.exp(1j * angles))))


def find_sync_rising_edge(t, x):
    """True buffer-wrap moment (n=0), located directly from the pulse
    data rather than assumed from the scope's trigger point.

    The threshold is taken from the median-to-maximum span, NOT from a
    percentile-based amplitude: the pulse occupies only a handful of
    samples in the record (40 us out of 16 ms is ~0.25%), so any
    percentile at or below the 99.5th sits underneath the pulse
    entirely and would yield a threshold down in the noise."""
    baseline = float(np.median(x))
    peak = float(np.max(x))
    span = peak - baseline

    if span < 0.5:
        raise ValueError(
            f"no clear sync pulse in CH2: peak is only {span*1e3:.1f} mV "
            f"above baseline. Check the probe is on PB0 and that this is "
            f"the CH2 file, not CH1.")

    thresh = baseline + 0.5 * span

    idx = np.where((x[:-1] < thresh) & (x[1:] >= thresh))[0]
    if len(idx) == 0:
        raise ValueError("no sync pulse rising edge found in CH2 data")

    i = int(idx[0])
    frac = (thresh - x[i]) / (x[i + 1] - x[i])
    t_sync = float(t[i] + frac * (t[i + 1] - t[i]))

    print(f"  CH2 baseline {baseline*1e3:.1f} mV, peak {peak*1e3:.1f} mV, "
          f"threshold {thresh*1e3:.1f} mV")
    print(f"  {len(idx)} rising edge(s) found in record")

    return t_sync


def extract_phase(t, x, fc, bw_hz, trim_cycles=10):
    fs = 1.0 / ((t[-1] - t[0]) / (len(t) - 1))
    x_ac = x - np.mean(x)
    z = hilbert(x_ac)
    zb = z * np.exp(-1j * 2 * np.pi * fc * t)
    b, a = butter(4, bw_hz / (fs / 2.0))
    zb = filtfilt(b, a, zb.real) + 1j * filtfilt(b, a, zb.imag)
    theta = np.unwrap(np.angle(zb))
    env = np.abs(zb)
    n = int(trim_cycles * fs / fc)
    sl = slice(n, -n)
    return t[sl], theta[sl], env[sl]


def bpsk_ideal_phase(t, t_sync, Rs, bits):
    k = np.floor((t - t_sync) * Rs).astype(int) % len(bits)
    return np.where(np.asarray(bits)[k] == 0, 0.0, np.pi)


def steady_state_mask(t, t_sync, Rs, guard_s=400e-6):
    phase_in_symbol = ((t - t_sync) * Rs) % 1.0
    edge_dist = np.minimum(phase_in_symbol, 1.0 - phase_in_symbol) / Rs
    return edge_dist > guard_s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ch1_csv")
    ap.add_argument("ch2_csv")
    ap.add_argument("--fc", type=float, required=True)
    ap.add_argument("--bw", type=float, default=BW_HZ)
    args = ap.parse_args()

    t1, x1 = load_tbs_csv(args.ch1_csv, channel=1)
    t2, x2 = load_tbs_csv(args.ch2_csv, channel=1)

    if len(t1) != len(t2):
        print("WARNING: CH1 and CH2 have different record lengths -- "
              "were they captured in the same acquisition?")

    t_sync = find_sync_rising_edge(t2, x2)
    print(f"Sync pulse found at t = {t_sync*1e3:.4f} ms")

    t_p, theta, env = extract_phase(t1, x1, args.fc, args.bw, TRIM_CYCLES)

    ideal = bpsk_ideal_phase(t_p, t_sync, RS, bits)
    mask = steady_state_mask(t_p, t_sync, RS, GUARD)

    err = wrap(theta - ideal)
    sym_idx = np.floor((t_p - t_sync) * RS).astype(int) % len(bits)
    state0 = mask & (np.asarray(bits)[sym_idx] == 0)
    state1 = mask & (np.asarray(bits)[sym_idx] == 1)

    if state0.sum() < 10 or state1.sum() < 10:
        print(f"WARNING: too few guard-excluded samples "
              f"(state0={state0.sum()}, state1={state1.sum()}) -- "
              f"check t_sync and Rs")

    mean0_deg = float(np.degrees(circular_mean(theta[state0])))
    mean1_deg = float(np.degrees(circular_mean(theta[state1])))
    separation = ((mean1_deg - mean0_deg + 180) % 360) - 180

    rms_deg = float(np.degrees(np.sqrt(np.mean(err[mask] ** 2))))
    peak_deg = float(np.degrees(np.max(np.abs(err[mask]))))
    env_ss = env[mask]
    ripple_ss_pct = float((env_ss.max() - env_ss.min()) / env_ss.mean() * 100)

    print(f"\nMean phase, state 0  : {mean0_deg:+.2f} deg  (ideal 0)")
    print(f"Mean phase, state 1  : {mean1_deg:+.2f} deg  (ideal 180)")
    print(f"Separation            : {separation:+.2f} deg  (ideal 180)")
    print(f"Steady-state RMS error: {rms_deg:.3f} deg")
    print(f"Peak error (ss only)  : {peak_deg:.3f} deg")
    print(f"Envelope ripple (ss)  : {ripple_ss_pct:.2f} %  (limit 5%)")
    print(f"Guard-excluded samples: state0={int(state0.sum())}, "
          f"state1={int(state1.sum())}")

    warn = []
    if not (150 <= abs(separation) <= 210):
        warn.append("separation is far from 180 deg -- check t_sync, the "
                     "bit pattern, or Rs")
    if ripple_ss_pct > 5:
        warn.append("steady-state envelope ripple exceeds 5%")

    print("\n=== Verdict ===")
    if warn:
        for w in warn:
            print(f"  WARNING: {w}")
    else:
        print("PASS -- looks like a good M3 capture.")


if __name__ == "__main__":
    main()
