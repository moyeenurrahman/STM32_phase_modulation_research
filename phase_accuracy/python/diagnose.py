"""
diagnose.py
============
Sanity-checks a TBS1000C CSV export (from a USB flash drive or from
OpenChoice) before it is trusted for analysis. Catches the failure modes
seen so far in this project:

  1. Wrong sample interval / record length -- wrong timebase, or a
     PC-software transfer that grabbed a live/partial buffer instead of
     the finished average.
  2. Garbled or noise-corrupted data -- transfer glitch, wrong channel
     selected, or the acquisition was read before averaging finished.

Works on the USB-drive CSV format and the OpenChoice CSV format alike:
the loader ignores any header rows (settings, model number, etc.) and
finds the numeric time/voltage columns automatically.

Usage:
    python diagnose.py ref_carrier_kp0.csv
    python diagnose.py m1_kp1p00_r1.csv --fc 5000 --interval 8e-6 --record 2000
    python diagnose.py ref_lpf_offset.csv --fc 5000 --channels 1,2

Exits 0 on PASS, 1 if any warning was raised. Do not proceed to analysis
on anything other than a clean PASS.
"""

import argparse
import sys
import numpy as np
from measure_fc_ref import load_tbs_csv        # reuse the existing loader


def sample_interval(t):
    return (t[-1] - t[0]) / (len(t) - 1)


def check_timebase(t, expected_interval, expected_record, tol=0.02):
    warnings = []
    dt = sample_interval(t)
    n = len(t)

    print(f"  record length        : {n}  (expect {expected_record})")
    print(f"  sample interval       : {dt*1e6:.4f} us  "
          f"(expect {expected_interval*1e6:.4f} us)")

    if abs(dt - expected_interval) / expected_interval > tol:
        warnings.append(
            f"sample interval off by {(dt/expected_interval-1)*100:+.1f}% "
            f"-- wrong timebase, or a partial/live acquisition was read")

    if abs(n - expected_record) > max(2, 0.01 * expected_record):
        warnings.append(
            f"record length is {n}, expected {expected_record} -- export "
            f"may be truncated or from the wrong source")

    d = np.diff(t)
    jitter_ppm = np.std(d) / np.mean(d) * 1e6
    print(f"  timebase grid jitter  : {jitter_ppm:.1f} ppm")
    if jitter_ppm > 500:
        warnings.append(
            "time column is irregular -- CSV may be corrupted or columns "
            "misaligned")

    if not np.all(d > 0):
        warnings.append(
            "time column is not strictly increasing -- rows are "
            "corrupted or out of order")

    return warnings


def check_channel(t, x, label, fc_expected):
    warnings = []
    xa = x - np.median(x)
    amp = (np.percentile(xa, 99.5) - np.percentile(xa, 0.5)) / 2.0
    dc = np.median(x)

    print(f"\n  --- {label} ---")
    print(f"  amplitude (peak)      : {amp:.4f}")
    print(f"  DC level              : {dc:.4f}")

    if amp < 1e-6:
        warnings.append(
            f"{label}: amplitude is essentially zero -- wrong channel, "
            f"no signal present, or the DAC/DMA is not running")
        return warnings

    # clipping / rail-pegged samples
    top = np.percentile(xa, 99.9)
    bot = np.percentile(xa, 0.1)
    flat_top = np.mean(xa > 0.995 * top)
    flat_bot = np.mean(xa < 0.995 * bot)
    print(f"  samples at extremes   : top {flat_top*100:.2f}%, "
          f"bottom {flat_bot*100:.2f}%")
    if flat_top > 0.02 or flat_bot > 0.02:
        warnings.append(f"{label}: possible clipping -- check vertical scale")

    # repeated / stuck values -- a classic sign of a corrupted transfer
    _, counts = np.unique(x, return_counts=True)
    max_run = int(counts.max())
    if max_run > 0.10 * len(x):
        warnings.append(
            f"{label}: {max_run} samples ({max_run/len(x)*100:.1f}%) share "
            f"one identical value -- likely a stuck channel or a "
            f"corrupted export, not real data")

    # implausible slew: compare the actual max |dx/dt| to what a clean
    # sine of this amplitude and fc could physically produce
    dt = sample_interval(t)
    dxdt = np.diff(x) / dt
    max_slew = float(np.max(np.abs(dxdt)))
    slew_limit = 2 * np.pi * fc_expected * amp * 3.0   # x3 safety margin
    print(f"  max |dx/dt|           : {max_slew:.3g}  "
          f"(physical limit ~{slew_limit:.3g})")
    if max_slew > slew_limit:
        warnings.append(
            f"{label}: sample-to-sample jumps exceed what a "
            f"{fc_expected:.0f} Hz sine of this amplitude can produce -- "
            f"data looks garbled, not filtered noise")

    # spectral concentration: how much power sits near fc vs everywhere else
    N = len(xa)
    w = np.blackman(N)
    X = np.abs(np.fft.rfft(xa * w)) ** 2
    f = np.fft.rfftfreq(N, dt)
    band = (f > fc_expected * 0.8) & (f < fc_expected * 1.2)
    frac_in_band = float(X[band].sum() / X.sum())
    print(f"  power near {fc_expected:.0f} Hz     : {frac_in_band*100:.1f}% "
          f"of total")
    if frac_in_band < 0.3:
        warnings.append(
            f"{label}: less than 30% of signal power is near "
            f"{fc_expected:.0f} Hz -- this looks like noise, not a "
            f"captured sine. Re-export or re-capture")

    return warnings


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--fc", type=float, default=5000.0,
                    help="expected dominant frequency, Hz")
    ap.add_argument("--interval", type=float, default=8e-6,
                    help="expected sample interval, s (8e-6 for 1 ms/div)")
    ap.add_argument("--record", type=int, default=2000,
                    help="expected record length, points")
    ap.add_argument("--channels", default="1",
                    help="comma-separated channel numbers to check, e.g. 1,2")
    args = ap.parse_args()

    channels = [int(c) for c in args.channels.split(",")]
    all_warnings = []

    print(f"=== {args.csv} ===")

    t_ref = None
    for ch in channels:
        t, x = load_tbs_csv(args.csv, channel=ch)
        if t_ref is None:
            print("\nTimebase:")
            all_warnings += check_timebase(t, args.interval, args.record)
            t_ref = t
        else:
            if len(t) != len(t_ref) or not np.allclose(t, t_ref, rtol=1e-6):
                all_warnings.append(
                    f"CH{ch} time column differs from CH1's -- columns may "
                    f"be misaligned in the CSV")
        all_warnings += check_channel(t, x, f"CH{ch}", args.fc)

    print("\n=== Verdict ===")
    if not all_warnings:
        print("PASS -- no problems detected. Safe to use for analysis.")
        sys.exit(0)
    else:
        for w in all_warnings:
            print(f"  WARNING: {w}")
        print(f"\n{len(all_warnings)} issue(s) found. Do not trust this "
              f"capture for analysis -- re-export or re-capture first.")
        sys.exit(1)


if __name__ == "__main__":
    main()
