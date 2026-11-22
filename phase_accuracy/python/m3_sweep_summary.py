"""
m3_sweep_summary.py
=====================
Aggregates multiple M3 (BPSK) capture repeats into the headline result:
mean phase per state, separation, steady-state RMS/peak error, and
envelope ripple -- averaged across repeats. Matches Table 10.4 (minus
the transition-time/overshoot columns, which come from the separate
zoomed capture, not from this whole-buffer analysis).

Reuses find_sync_rising_edge, extract_phase, bpsk_ideal_phase, and
steady_state_mask from check_m3_point.py.

Usage:
    Create a manifest text file, one repeat per line:
        "ch1_file  ch2_file"

        m3_ch1_r1.csv  m3_ch2_r1.csv
        m3_ch1_r2.csv  m3_ch2_r2.csv
        m3_ch1_r3.csv  m3_ch2_r3.csv

    python m3_sweep_summary.py m3_manifest.txt --fc 5000.12
"""

import argparse
import numpy as np
from check_m3_point import (
    find_sync_rising_edge, extract_phase, bpsk_ideal_phase,
    steady_state_mask, wrap, circular_mean, RS, GUARD, BW_HZ,
    TRIM_CYCLES, bits,
)
from measure_fc_ref import load_tbs_csv


def analyze_one(ch1_csv, ch2_csv, fc, bw):
    t1, x1 = load_tbs_csv(ch1_csv, channel=1)
    t2, x2 = load_tbs_csv(ch2_csv, channel=1)

    t_sync = find_sync_rising_edge(t2, x2)
    t_p, theta, env = extract_phase(t1, x1, fc, bw, TRIM_CYCLES)

    ideal = bpsk_ideal_phase(t_p, t_sync, RS, bits)
    mask = steady_state_mask(t_p, t_sync, RS, GUARD)
    err = wrap(theta - ideal)

    sym_idx = np.floor((t_p - t_sync) * RS).astype(int) % len(bits)
    state0 = mask & (np.asarray(bits)[sym_idx] == 0)
    state1 = mask & (np.asarray(bits)[sym_idx] == 1)

    mean0 = float(np.degrees(circular_mean(theta[state0])))
    mean1 = float(np.degrees(circular_mean(theta[state1])))
    separation = ((mean1 - mean0 + 180) % 360) - 180

    rms = float(np.degrees(np.sqrt(np.mean(err[mask] ** 2))))
    peak = float(np.degrees(np.max(np.abs(err[mask]))))
    env_ss = env[mask]
    ripple = float((env_ss.max() - env_ss.min()) / env_ss.mean() * 100)

    return {
        "t_sync": t_sync, "mean0": mean0, "mean1": mean1,
        "separation": separation, "rms": rms, "peak": peak,
        "ripple": ripple, "n0": int(state0.sum()), "n1": int(state1.sum()),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("manifest")
    ap.add_argument("--fc", type=float, required=True)
    ap.add_argument("--bw", type=float, default=BW_HZ)
    args = ap.parse_args()

    results = []
    with open(args.manifest) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) != 2:
                print(f"WARNING: skipping malformed line: {line}")
                continue
            ch1, ch2 = parts
            r = analyze_one(ch1, ch2, args.fc, args.bw)
            results.append(r)
            flag = "  *** RIPPLE > 5% ***" if r["ripple"] > 5 else ""
            print(f"[{ch1}]")
            print(f"  separation={r['separation']:+.2f} deg  "
                  f"RMS={r['rms']:.2f} deg  peak={r['peak']:.2f} deg  "
                  f"ripple={r['ripple']:.2f}%{flag}  "
                  f"(n0={r['n0']}, n1={r['n1']})")

    if not results:
        print("No valid rows found in manifest.")
        return

    def agg(key):
        vals = np.array([r[key] for r in results])
        return float(vals.mean()), float(vals.std())

    def circular_agg_deg(vals_deg):
        angles = np.radians(np.array(vals_deg))
        C = np.mean(np.cos(angles))
        S = np.mean(np.sin(angles))
        R = np.hypot(C, S)
        mean_deg = float(np.degrees(np.arctan2(S, C)))
        circ_std_deg = (float(np.degrees(np.sqrt(max(-2 * np.log(R), 0))))
                        if R > 1e-9 else float("nan"))
        return mean_deg, circ_std_deg

    # separation, mean0, mean1 are ANGLES -- must use circular statistics,
    # not plain mean/std, or values straddling +/-180 average to garbage
    sep_m, sep_s = circular_agg_deg([r["separation"] for r in results])
    m0_m, m0_s = circular_agg_deg([r["mean0"] for r in results])
    m1_m, m1_s = circular_agg_deg([r["mean1"] for r in results])
    # RMS, peak, ripple are magnitudes -- plain mean/std is correct
    rms_m, rms_s = agg("rms")
    peak_m, peak_s = agg("peak")
    ripple_m, ripple_s = agg("ripple")

    print(f"\n=== M3 summary ({len(results)} repeats) ===")
    print(f"Mean phase, state 0    : {m0_m:+.2f} +/- {m0_s:.2f} deg")
    print(f"Mean phase, state 1    : {m1_m:+.2f} +/- {m1_s:.2f} deg")
    print(f"Separation              : {sep_m:+.2f} +/- {sep_s:.2f} deg  "
          f"(ideal 180)")
    print(f"Steady-state RMS error  : {rms_m:.2f} +/- {rms_s:.2f} deg")
    print(f"Peak error (ss only)    : {peak_m:.2f} +/- {peak_s:.2f} deg")
    print(f"Envelope ripple (ss)    : {ripple_m:.2f} +/- {ripple_s:.2f} %")
    print(f"\nGuard interval used: +/-{GUARD*1e6:.0f} us")
    print(f"Post-detection bandwidth: {args.bw:.0f} Hz")

    if abs(abs(sep_m) - 180) > 5:
        print("\nWARNING: mean separation is more than 5 deg from 180 -- "
              "check t_sync detection and the bit pattern across repeats.")
    if ripple_m > 5:
        print("\nWARNING: mean steady-state ripple exceeds 5%.")


if __name__ == "__main__":
    main()
