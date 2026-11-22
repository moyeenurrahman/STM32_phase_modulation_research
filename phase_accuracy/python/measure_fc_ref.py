"""
measure_fc_ref.py
=================
Analyse the R1 reference capture (ref_carrier_kp0.csv) and produce the
measured carrier frequency f_c_hat used by every phase measurement.

Usage:
    python measure_fc_ref.py ref_carrier_kp0.csv
    python measure_fc_ref.py ref_carrier_kp0.csv --channel 1

Outputs:
    - printed health checks and three independent frequency estimates
    - fc_ref.json, read automatically by the phase analysis scripts

Three estimators are used. They should agree. Disagreement means a
problem with the capture, not with the modulator.

    1. Zero-crossing regression   -- primary, sub-ppm
    2. Hilbert phase slope        -- cross-check, sub-ppm
    3. FFT peak, parabolic fit    -- coarse sanity check, ~100 ppm

NOTE: if analyse_freq.py exposes an equivalent estimator, run it too and
confirm agreement to better than 0.05 ppm.
"""

import sys
import json
import argparse
import numpy as np
from scipy.signal import hilbert

FC_NOMINAL = 5000.0
FS_DAC = 100000.0          # DAC update rate
CYCLES_PER_BUFFER = 50     # 1000 samples at 100 kSa/s, fc = fs/20


# ----------------------------------------------------------------------
# CSV loading
# ----------------------------------------------------------------------
def load_tbs_csv(path, channel=1):
    """
    Robust loader for TBS1000C CSV exports.

    Handles the metadata block in the leading columns and returns
    (time, voltage) for the requested channel. channel=1 is the first
    waveform column, channel=2 the second.
    """
    rows = []
    with open(path, "r", errors="ignore") as f:
        for line in f:
            parts = [p.strip() for p in line.split(",")]
            vals = []
            for p in parts:
                try:
                    vals.append(float(p))
                except ValueError:
                    vals.append(np.nan)
            if np.count_nonzero(~np.isnan(vals)) >= 2:
                rows.append(vals)

    width = max(len(r) for r in rows)
    data = np.full((len(rows), width), np.nan)
    for i, r in enumerate(rows):
        data[i, :len(r)] = r

    # keep columns that are numeric almost everywhere
    good = np.mean(~np.isnan(data), axis=0) > 0.95
    data = data[:, good]
    data = data[~np.isnan(data).any(axis=1)]

    if data.shape[1] < 2:
        raise ValueError("Could not find time and waveform columns")

    # time column: mostly monotonically increasing, with a sensible
    # positive median step. Selection relies only on monotonicity, not
    # on sample-to-sample step uniformity: real scope CSV exports can
    # show per-sample rounding jitter in the printed time values at
    # fine timebases (large absolute time values relative to a small
    # step), even though the actual acquisition clock is genuinely
    # uniform. This doesn't affect correctness downstream, since every
    # script computes the sample interval via endpoint division over
    # the whole record, not from individual consecutive differences.
    tcol = None
    best_score = -1.0
    for c in range(data.shape[1]):
        d = np.diff(data[:, c])
        if len(d) == 0:
            continue
        frac_positive = np.mean(d > 0)
        med = np.median(d)
        if med <= 0:
            continue
        if frac_positive > 0.99 and frac_positive > best_score:
            best_score = frac_positive
            tcol = c
    if tcol is None:
        raise ValueError("Could not identify the time column")

    vcols = [c for c in range(data.shape[1]) if c != tcol]
    if channel > len(vcols):
        raise ValueError(f"Only {len(vcols)} waveform column(s) in this file")

    return data[:, tcol], data[:, vcols[channel - 1]]


def sample_interval(t):
    """Endpoint division -- robust to limited precision in the CSV."""
    return (t[-1] - t[0]) / (len(t) - 1)


# ----------------------------------------------------------------------
# Health checks
# ----------------------------------------------------------------------
def health_checks(t, x):
    dt = sample_interval(t)
    fs = 1.0 / dt
    dur = t[-1] - t[0]

    print("--- Capture health ---")
    print(f"  points              : {len(t)}")
    print(f"  sample interval     : {dt*1e6:.3f} us   (target 8.000 us)")
    print(f"  sample rate         : {fs/1e3:.3f} kSa/s")
    print(f"  record duration     : {dur*1e3:.3f} ms")
    print(f"  carrier cycles      : {dur*FC_NOMINAL:.1f}")
    print(f"  samples per cycle   : {fs/FC_NOMINAL:.2f}")

    # timing jitter in the CSV itself
    d = np.diff(t)
    print(f"  timebase step spread: {np.std(d)/np.mean(d)*1e6:.2f} ppm")

    # amplitude and clipping
    xa = x - np.median(x)
    amp = (np.percentile(xa, 99.5) - np.percentile(xa, 0.5)) / 2.0
    top = np.percentile(xa, 99.9)
    bot = np.percentile(xa, 0.1)
    flat_top = np.mean(xa > 0.995 * top)
    flat_bot = np.mean(xa < 0.995 * bot)
    print(f"  amplitude           : {amp:.4f} V (peak)")
    print(f"  samples at extremes : top {flat_top*100:.2f} %, "
          f"bottom {flat_bot*100:.2f} %")

    warn = []
    if abs(dt - 8e-6) / 8e-6 > 0.02:
        warn.append("sample interval is not 8 us -- wrong timebase")
    if fs / FC_NOMINAL < 20:
        warn.append("fewer than 20 samples per carrier cycle")
    if dur * FC_NOMINAL < 40:
        warn.append("fewer than 40 carrier cycles in the record")
    if flat_top > 0.02 or flat_bot > 0.02:
        warn.append("waveform may be clipped -- check vertical scale")
    return warn


# ----------------------------------------------------------------------
# Estimator 1: zero-crossing regression with Schmitt hysteresis
# ----------------------------------------------------------------------
def zc_frequency(t, x, hyst_frac=0.15, trim_frac=0.01, passes=3):
    xa = x - np.median(x)
    amp = (np.percentile(xa, 99.5) - np.percentile(xa, 0.5)) / 2.0
    h = hyst_frac * amp

    # rising crossings, armed only after the signal has gone below -h
    t0 = []
    armed = False
    for i in range(len(xa) - 1):
        if xa[i] < -h:
            armed = True
        if armed and xa[i] < 0.0 <= xa[i + 1]:
            frac = -xa[i] / (xa[i + 1] - xa[i])
            t0.append(t[i] + frac * (t[i + 1] - t[i]))
            armed = False
    t0 = np.asarray(t0)

    if len(t0) < 10:
        raise ValueError(f"only {len(t0)} zero crossings found")

    # trim the record edges
    ntrim = max(1, int(np.ceil(trim_frac * len(t0))))
    t0 = t0[ntrim:len(t0) - ntrim]

    # time-based cycle indexing with iterative refinement
    f = 1.0 / np.mean(np.diff(t0))
    for _ in range(passes):
        k = np.round((t0 - t0[0]) * f)
        A = np.vstack([k, np.ones_like(k)]).T
        coef, *_ = np.linalg.lstsq(A, t0, rcond=None)
        f = 1.0 / coef[0]

    resid = t0 - (A @ coef)
    jitter_rms = float(np.std(resid))

    # standard error of the slope -> frequency uncertainty
    n = len(k)
    s2 = np.sum(resid ** 2) / (n - 2)
    se_slope = np.sqrt(s2 / np.sum((k - k.mean()) ** 2))
    sigma_f = f ** 2 * se_slope

    return {
        "f": float(f),
        "n_crossings": int(n),
        "jitter_rms_ns": jitter_rms * 1e9,
        "sigma_f_hz": float(sigma_f),
        "sigma_f_ppm": float(sigma_f / f * 1e6),
        "max_resid_ns": float(np.max(np.abs(resid)) * 1e9),
    }


# ----------------------------------------------------------------------
# Estimator 2: Hilbert phase slope
# ----------------------------------------------------------------------
def phase_slope_frequency(t, x, trim_cycles=10, f_guess=FC_NOMINAL):
    fs = 1.0 / sample_interval(t)
    xa = x - np.mean(x)
    z = hilbert(xa)
    ph = np.unwrap(np.angle(z))

    n = int(trim_cycles * fs / f_guess)
    sl = slice(n, len(t) - n)
    tt, pp = t[sl], ph[sl]

    A = np.vstack([tt, np.ones_like(tt)]).T
    coef, *_ = np.linalg.lstsq(A, pp, rcond=None)
    return float(coef[0] / (2.0 * np.pi))


# ----------------------------------------------------------------------
# Estimator 3: FFT peak with parabolic interpolation
# ----------------------------------------------------------------------
def fft_frequency(t, x):
    dt = sample_interval(t)
    xa = x - np.mean(x)
    N = len(xa)
    w = np.blackman(N)
    X = np.abs(np.fft.rfft(xa * w))
    f = np.fft.rfftfreq(N, dt)

    k = int(np.argmax(X[1:]) + 1)
    a, b, c = np.log(X[k - 1]), np.log(X[k]), np.log(X[k + 1])
    delta = 0.5 * (a - c) / (a - 2 * b + c)
    return float((k + delta) * (f[1] - f[0]))


# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--channel", type=int, default=1,
                    help="1 = first waveform column (LPF output)")
    ap.add_argument("--out", default="fc_ref.json")
    args = ap.parse_args()

    t, x = load_tbs_csv(args.csv, args.channel)
    warn = health_checks(t, x)

    print("\n--- Frequency estimates ---")
    zc = zc_frequency(t, x)
    f_zc = zc["f"]
    print(f"  1. zero crossing    : {f_zc:.6f} Hz "
          f"({(f_zc/FC_NOMINAL - 1)*1e6:+.3f} ppm vs nominal)")
    print(f"       crossings used : {zc['n_crossings']}")
    print(f"       fit jitter RMS : {zc['jitter_rms_ns']:.1f} ns")
    print(f"       worst residual : {zc['max_resid_ns']:.1f} ns")
    print(f"       1-sigma        : {zc['sigma_f_ppm']:.4f} ppm")

    f_ps = phase_slope_frequency(t, x)
    print(f"  2. phase slope      : {f_ps:.6f} Hz "
          f"({(f_ps/f_zc - 1)*1e6:+.3f} ppm vs estimator 1)")

    f_ft = fft_frequency(t, x)
    print(f"  3. FFT parabolic    : {f_ft:.3f} Hz "
          f"({(f_ft/f_zc - 1)*1e6:+.0f} ppm vs estimator 1)")

    # implied DAC sample rate, since fc = fs/20 by construction
    fs_implied = f_zc * (1000.0 / CYCLES_PER_BUFFER)

    print("\n--- Consistency ---")
    print(f"  implied DAC fs      : {fs_implied:.3f} Sa/s "
          f"({(fs_implied/FS_DAC - 1)*1e6:+.2f} ppm vs 100 kSa/s)")
    print("  Compare this ppm figure with the frequency-accuracy metric.")
    print("  They come from the same MCU clock and must agree, within the")
    print("  scope timebase error common to both.")

    if abs(f_ps - f_zc) / f_zc * 1e6 > 5.0:
        warn.append("estimators 1 and 2 differ by more than 5 ppm")
    if abs(f_ft - f_zc) / f_zc * 1e6 > 500.0:
        warn.append("FFT estimate is far off -- wrong column or bad capture")
    if zc["max_resid_ns"] > 20 * zc["jitter_rms_ns"]:
        warn.append("an outlier crossing was found -- possible glitch")

    print("\n--- Verdict ---")
    if warn:
        for w in warn:
            print(f"  WARNING: {w}")
    else:
        print("  All checks passed.")

    result = {
        "source_file": args.csv,
        "fc_hat_hz": f_zc,
        "sigma_ppm": zc["sigma_f_ppm"],
        "offset_ppm_vs_nominal": (f_zc / FC_NOMINAL - 1) * 1e6,
        "f_phase_slope_hz": f_ps,
        "f_fft_hz": f_ft,
        "implied_fs_hz": fs_implied,
        "sample_interval_s": sample_interval(t),
        "n_points": int(len(t)),
        "warnings": warn,
    }
    with open(args.out, "w") as f:
        json.dump(result, f, indent=2)

    print(f"\n  f_c_hat = {f_zc:.6f} Hz  ->  written to {args.out}")
    print("  Every phase script must read this file. Do not retype it.")


if __name__ == "__main__":
    main()
