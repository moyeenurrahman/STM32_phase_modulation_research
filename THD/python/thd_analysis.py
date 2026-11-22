#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
 THD ANALYSIS  --  BSc thesis, baseband phase modulation with STM32
===============================================================================

 Reads Tektronix TBS1000C CSV captures, computes total harmonic distortion
 (THD_F, harmonics 2..7, IEEE Std 1241 definition) and produces the tables
 and figures needed for the THD chapter.

 Written for 2000-point oscilloscope records, NOT for the 262144-point
 records assumed in the original methodology document.

 Usage
 -----
   1. Put every CSV in one folder.
   2. Set DATA_DIR below.
   3. Check that the file patterns in CASES match your filenames.
   4. python thd_analysis.py

 Output
 ------
   - printed tables for Tables 1, 2 and 3 of the THD document
   - figures/*.pdf and figures/*.png at 300 dpi
   - results_thd.csv  with every number, for pasting into the thesis

===============================================================================
"""

import os
import glob
import math
import hashlib
import numpy as np
import matplotlib
matplotlib.use("Agg")           # no display needed
import matplotlib.pyplot as plt
from scipy.signal.windows import blackmanharris

# =============================================================================
# 1.  CONFIGURATION  --  edit this section only
# =============================================================================

DATA_DIR = "."                  # folder holding the CSV files
FIG_DIR  = "figures"            # output folder for plots

N_HARMONICS   = 7               # count harmonics 2..7
LOBE_BINS     = 4               # Blackman-Harris main lobe half width, in bins
FLOOR_MARGIN_DB = 6.0           # harmonic must beat the floor by this much

LPF_FCUT      = 20100.0         # measured Sallen-Key cutoff, Hz
LPF_ORDER     = 2               # 2nd order Butterworth

VDDA          = 3.3             # for converting DAC codes to volts
DAC_FULLSCALE = 4095

# -----------------------------------------------------------------------------
# Case table.
#
#   pattern  : glob pattern that matches the CSV files for this case
#   f0       : nominal fundamental, Hz
#   codes    : DAC amplitude in codes (None for floor captures)
#   point    : "A", "B" or None
#   floor    : name of the case to use as the external noise floor reference
#   lpf_corr : True only where the filter response must be de-embedded
# -----------------------------------------------------------------------------

CASES = {
    # ---- noise floors (analysed first) --------------------------------------
    "T9a": dict(pattern="T9a/*.csv", f0=1000.0, codes=None, point="B",
                floor=None, lpf_corr=False, kind="floor",
                note="system floor, 500 mV/div, 4 ms/div"),

    "T9b": dict(pattern="T9b/*.csv", f0=1000.0, codes=None, point="B",
                floor=None, lpf_corr=False, kind="floor",
                note="system floor, 200 mV/div, 4 ms/div"),

    "T9c": dict(pattern="T9c/*.csv", f0=5000.0, codes=None, point="B",
                floor=None, lpf_corr=False, kind="floor",
                note="system floor, 500 mV/div, 400 us/div, from SNR N1"),

    # ---- instrument distortion check ---------------------------------------
    "FG":  dict(pattern="FG/*.csv", f0=1000.0, codes=None, point="-",
                floor="T9a", lpf_corr=False, kind="signal",
                note="function generator direct into scope, 1.93 Vpp"),

    # ---- primary results ----------------------------------------------------
    "T1":  dict(pattern="T1/*.csv", f0=1000.0, codes=1200, point="A",
                floor="T9a", lpf_corr=False, kind="signal",
                note="raw DAC linearity"),

    "T2":  dict(pattern="T2/*.csv", f0=1000.0, codes=1200, point="B",
                floor="T9a", lpf_corr=False, kind="signal",
                note="HEADLINE RESULT"),

    # ---- amplitude sweep, Point B ------------------------------------------
    "T5B": dict(pattern="T5B/*.csv", f0=1000.0, codes=500, point="B",
                floor="T9b", lpf_corr=False, kind="signal", note="sweep"),

    "T6B": dict(pattern="T6B/*.csv", f0=1000.0, codes=800, point="B",
                floor="T9b", lpf_corr=False, kind="signal", note="sweep"),

    "T7":  dict(pattern="T7/*.csv", f0=1000.0, codes=1000, point="B",
                floor="T9b", lpf_corr=False, kind="signal",
                note="sweep - CHECK Vpp, 200 mV/div may have clipped"),

    "T8B": dict(pattern="T8B/*.csv", f0=1000.0, codes=1500, point="B",
                floor="T9a", lpf_corr=False, kind="signal",
                note="flat tops, op-amp input handover, LOWER BOUND"),

    # ---- amplitude sweep, Point A ------------------------------------------
    "T5A": dict(pattern="T5A/*.csv", f0=1000.0, codes=500, point="A",
                floor="T9b", lpf_corr=False, kind="signal", note="sweep"),

    "T6A": dict(pattern="T6A/*.csv", f0=1000.0, codes=800, point="A",
                floor="T9b", lpf_corr=False, kind="signal", note="sweep"),

    "T8A": dict(pattern="T8A/*.csv", f0=1000.0, codes=1500, point="A",
                floor="T9a", lpf_corr=False, kind="signal",
                note="clean at Point A - proves the op-amp is the cause"),

    # ---- 5 kHz operating point (reused from SNR N2 / N3) --------------------
    "T3":  dict(pattern="T3/*.csv", f0=5000.0, codes=1200, point="A",
                floor="T9c", lpf_corr=False, kind="signal",
                note="operating point, raw, reused from SNR N2"),

    "T4":  dict(pattern="T4/*.csv", f0=5000.0, codes=1200, point="B",
                floor="T9c", lpf_corr=True, kind="signal",
                note="operating point, filter corrected, reused from SNR N3"),
}

# Which cases form the amplitude sweep plot
SWEEP_B = ["T5B", "T6B", "T7", "T2", "T8B"]
SWEEP_A = ["T5A", "T6A", "T1", "T8A"]

# Order in which cases are analysed and reported
RUN_ORDER = ["T9a", "T9b", "T9c", "FG",
             "T2", "T1",
             "T5B", "T6B", "T7", "T8B",
             "T5A", "T6A", "T8A",
             "T3", "T4"]

# For the matched-run comparison against the 1 kHz cases
MATCHED_SUBSET_N = 8


# =============================================================================
# 2.  CSV LOADING
# =============================================================================

def load_scope_csv(path):
    """
    Read a Tektronix CSV and return (time, volts) as float arrays.

    Handles both the pendrive format and the OpenChoice format by ignoring
    the header and looking for the numeric block.  The time column is
    identified as the column with a near-constant positive step.
    """
    rows = []
    header = {}

    with open(path, "r", errors="ignore") as fh:
        for line in fh:
            parts = [p.strip() for p in line.rstrip("\n").split(",")]

            # keep header keywords such as "Sample Interval"
            if len(parts) >= 2 and parts[0] and not _is_float(parts[0]):
                if _is_float(parts[1]):
                    header[parts[0].lower()] = float(parts[1])

            vals = [float(p) if _is_float(p) else np.nan for p in parts]
            if sum(not math.isnan(v) for v in vals) >= 2:
                rows.append(vals)

    if not rows:
        raise ValueError(f"No numeric data found in {path}")

    width = max(len(r) for r in rows)
    arr = np.full((len(rows), width), np.nan)
    for i, r in enumerate(rows):
        arr[i, :len(r)] = r

    # keep only columns with no missing values
    good = [c for c in range(width) if np.all(np.isfinite(arr[:, c]))]
    if len(good) < 1:
        raise ValueError(f"Could not find clean columns in {path}")

    # ---- identify the time column ---------------------------------------
    # A time column rises strictly from start to finish.  A voltage column
    # does not.  That single fact separates them reliably, and it matters:
    # on a capture with no signal the voltage column has a very small range,
    # so picking "the column with the largest range" would grab the time
    # column instead and analyse the timebase as if it were the waveform.
    strictly_rising = []
    for c in good:
        d = np.diff(arr[:, c])
        if len(d) and np.all(d > 0):
            strictly_rising.append(c)

    t_col = None
    for c in strictly_rising:
        d = np.diff(arr[:, c])
        if np.std(d) <= 1e-3 * abs(np.mean(d)):      # near constant step
            t_col = c
            break

    # the voltage column must NOT be strictly rising
    cands = [c for c in good if c not in strictly_rising]
    if not cands:
        cands = [c for c in good if c != t_col]
    if not cands:
        raise ValueError(f"No voltage column in {path}")
    v_col = max(cands, key=lambda c: np.ptp(arr[:, c]))

    v = arr[:, v_col]

    if t_col is not None:
        t = arr[:, t_col]
    else:
        dt = header.get("sample interval")
        if dt is None:
            raise ValueError(f"No time column and no Sample Interval in {path}")
        t = np.arange(len(v)) * dt

    return t, v


def unique_files(files):
    """
    Drop byte-identical duplicates.

    Saving a capture twice without pressing Single first produces two files
    with the same contents.  They are not two measurements: keeping both
    makes the run-to-run spread look smaller than it really is.
    """
    seen = {}
    out = []
    dropped = []
    for p in files:
        h = hashlib.md5(open(p, "rb").read()).hexdigest()
        if h in seen:
            dropped.append(os.path.basename(p))
        else:
            seen[h] = p
            out.append(p)
    return out, dropped


def _is_float(s):
    try:
        float(s)
        return True
    except (TypeError, ValueError):
        return False


# =============================================================================
# 3.  SPECTRUM AND HARMONIC EXTRACTION
# =============================================================================

def spectrum(x, dt):
    """
    One-sided power spectrum using a Blackman-Harris window.

    Blackman-Harris is used rather than Hann because the STM32 clock and the
    oscilloscope clock are independent.  Hann side lobes near -31 dB would
    bury harmonics sitting 50 to 60 dB below the fundamental.

    The window normalisation cancels when powers are ratioed, so no window
    correction factor is applied.
    """
    x = np.asarray(x, dtype=float)
    x = x - np.mean(x)                       # remove DC before the FFT
    n = len(x)
    w = blackmanharris(n)
    X = np.fft.rfft(x * w)
    f = np.fft.rfftfreq(n, d=dt)
    P = np.abs(X) ** 2
    return f, P


def lobe_power(P, df, f_target):
    """
    Sum the power across the whole window main lobe.

    A single FFT bin holds only part of a component's energy because the
    Blackman-Harris main lobe spans about 8 bins.  Taking the maximum of one
    bin underestimates every harmonic and makes THD look better than it is.
    """
    k = int(round(f_target / df))
    lo = max(k - LOBE_BINS, 0)
    hi = min(k + LOBE_BINS + 1, len(P))
    if lo >= hi:
        return 0.0
    return float(np.sum(P[lo:hi]))


def local_floor(P, df, f_target, f0):
    """
    Estimate the noise floor beside a peak, scaled to the same number of bins
    as the main lobe.

    The guard and width are derived from the harmonic spacing, because a
    2000-point record gives only about 20 to 40 bins between harmonics.  A
    fixed width of 60 bins, as in the original document, would reach into the
    neighbouring harmonic and inflate the floor.
    """
    bins_between = f0 / df
    guard = LOBE_BINS + 2
    width = int(max(4, bins_between / 2.0 - guard - 1))

    k = int(round(f_target / df))
    left = P[max(k - guard - width, 0): max(k - guard, 1)]
    right = P[min(k + guard, len(P) - 1): min(k + guard + width, len(P))]
    ref = np.concatenate([left, right])
    if len(ref) == 0:
        return 0.0
    return float(np.median(ref) * (2 * LOBE_BINS + 1))


def butter_gain(freq, fcut=LPF_FCUT, order=LPF_ORDER):
    """Magnitude response of the Butterworth low-pass filter."""
    return 1.0 / math.sqrt(1.0 + (freq / fcut) ** (2 * order))


def analyse_run(path, f0_nom, external_floor=None, lpf_corr=False):
    """
    Analyse one CSV file.  Returns a dict of per-harmonic results.

    external_floor : optional dict {harmonic_number: floor_power} taken from
                     a T9 capture at the same volts/div and time/div.
    """
    t, v = load_scope_csv(path)
    dt = float(np.mean(np.diff(t)))
    fs = 1.0 / dt
    n = len(v)

    f, P = spectrum(v, dt)
    df = f[1] - f[0]

    # locate the true fundamental, which will not be exactly f0_nom because
    # the STM32 and scope clocks are independent
    band = np.abs(f - f0_nom) <= max(50.0, 3 * df)
    if not np.any(band):
        raise ValueError(f"No bins near {f0_nom} Hz in {path}")
    k0 = np.where(band)[0][int(np.argmax(P[band]))]
    f0 = f[k0]

    P1 = lobe_power(P, df, f0)
    g1 = butter_gain(f0) if lpf_corr else 1.0
    P1c = P1 / g1 ** 2

    harmonics = {}
    for h in range(2, N_HARMONICS + 1):
        fh = h * f0
        if fh >= fs / 2.0:
            harmonics[h] = dict(freq=fh, ratio=0.0, status="above Nyquist",
                                measured=False)
            continue

        Ph = lobe_power(P, df, fh)
        fl = local_floor(P, df, fh, f0)

        # use the higher of the local estimate and the external T9 floor
        if external_floor is not None and h in external_floor:
            fl = max(fl, external_floor[h])

        threshold = fl * (10 ** (FLOOR_MARGIN_DB / 10.0))
        measured = Ph > threshold

        Ph_net = max(Ph - fl, 0.0)
        use = Ph_net if measured else fl        # floor acts as an upper bound

        if lpf_corr:
            gh = butter_gain(fh)
            use = use / gh ** 2

        harmonics[h] = dict(freq=fh,
                            ratio=use / P1c if P1c > 0 else 0.0,
                            floor_ratio=(fl / (butter_gain(fh) ** 2 if lpf_corr else 1.0)) / P1c if P1c > 0 else 0.0,
                            status="measured" if measured else "at/below floor",
                            measured=measured)

    total = sum(v_["ratio"] for v_ in harmonics.values())
    thd_f = math.sqrt(total) if total > 0 else 0.0

    return dict(path=os.path.basename(path), f0=f0, fs=fs, n=n, df=df,
                P1=P1c, harmonics=harmonics, thd_f=thd_f,
                f_axis=f, P_axis=P, vpp=float(np.ptp(v)),
                crest_pos=crest_factors(v)[0],
                crest_neg=crest_factors(v)[1])


def crest_factors(v):
    """
    Return the positive and negative crest factors of the record.

    Crest factor is peak divided by RMS.  For an undistorted sine it is
    exactly sqrt(2) = 1.414, whatever the amplitude, sample rate or record
    length.  Flattening a peak removes area from the extreme of the
    waveform, so the crest factor on that side falls.

    This replaces the earlier test that counted samples near the peak.  That
    test failed on real captures: an 8-bit oscilloscope filling half the
    screen resolves the signal in steps of about 16 mV, which is larger than
    the 1 % slices the test relied on, so the slices were often empty and
    the ratio came out as infinity.  Crest factor uses every sample and is
    immune to the code step.

    Guide values
    ------------
      1.414          pure sine
      below 1.36     that peak is flattened
      pos < neg      asymmetric, e.g. an op-amp input handover
      above 1.45     noise or spikes dominate, not a clean sine
    """
    v = np.asarray(v, dtype=float)
    v = v - np.mean(v)
    r = float(np.std(v))
    if r <= 0:
        return float("nan"), float("nan")
    return float(np.max(v) / r), float(-np.min(v) / r)


def expected_vpp(codes):
    """Theoretical peak-to-peak output for a given DAC amplitude in codes."""
    if codes is None:
        return None
    return 2.0 * codes * VDDA / DAC_FULLSCALE


def analyse_floor_file(path, f0_nom):
    """
    Analyse a T9 capture and return the floor power at each harmonic
    frequency, in the same lobe-summed units used for the signal files.
    """
    t, v = load_scope_csv(path)
    dt = float(np.mean(np.diff(t)))
    f, P = spectrum(v, dt)
    df = f[1] - f[0]

    out = {}
    for h in range(2, N_HARMONICS + 1):
        fh = h * f0_nom
        if fh >= (1.0 / dt) / 2.0:
            continue
        out[h] = lobe_power(P, df, fh)
    return out, dict(f_axis=f, P_axis=P, path=os.path.basename(path))


# =============================================================================
# 4.  CASE-LEVEL AGGREGATION
# =============================================================================

def aggregate(runs):
    """
    Combine several runs of the same case.

    The per-harmonic POWER RATIOS are averaged, not the raw spectra.  The
    measured fundamental sits in a slightly different bin in every run, so
    averaging spectra directly would smear the harmonics and understate THD.
    Each run is analysed on its own detected fundamental first, which removes
    that problem.
    """
    if not runs:
        return None

    hs = range(2, N_HARMONICS + 1)
    mean_ratio = {}
    n_measured = {}
    for h in hs:
        vals = [r["harmonics"][h]["ratio"] for r in runs if h in r["harmonics"]]
        mean_ratio[h] = float(np.mean(vals)) if vals else 0.0
        n_measured[h] = sum(1 for r in runs
                            if h in r["harmonics"] and r["harmonics"][h]["measured"])

    total = sum(mean_ratio.values())
    thd_f = math.sqrt(total) if total > 0 else 0.0

    per_run = [r["thd_f"] for r in runs]
    thd_mean = float(np.mean(per_run))
    thd_std = float(np.std(per_run, ddof=1)) if len(per_run) > 1 else 0.0

    # spread in dB, useful for the uncertainty statement
    if thd_mean > 0:
        db_vals = [20 * math.log10(x) for x in per_run if x > 0]
        db_std = float(np.std(db_vals, ddof=1)) if len(db_vals) > 1 else 0.0
    else:
        db_std = 0.0

    dom = max(mean_ratio, key=lambda h: mean_ratio[h]) if total > 0 else None

    return dict(n_runs=len(runs),
                f0=float(np.mean([r["f0"] for r in runs])),
                fs=runs[0]["fs"], n=runs[0]["n"], df=runs[0]["df"],
                vpp=float(np.mean([r["vpp"] for r in runs])),
                crest_pos=float(np.median([r.get("crest_pos", np.nan) for r in runs])),
                crest_neg=float(np.median([r.get("crest_neg", np.nan) for r in runs])),
                ratio=mean_ratio, n_measured=n_measured,
                thd_f=thd_f, thd_mean=thd_mean, thd_std=thd_std,
                thd_db_std=db_std, dominant=dom, runs=runs)


def db(x):
    return 20 * math.log10(x) if x > 0 else float("nan")


def dbc(ratio):
    return 10 * math.log10(ratio) if ratio > 0 else float("nan")


# =============================================================================
# 5.  MAIN
# =============================================================================

def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    results = {}
    floors = {}
    floor_spectra = {}

    print("=" * 74)
    print(" THD ANALYSIS")
    print("=" * 74)

    # ---- pass 1: noise floors ------------------------------------------------
    for name in RUN_ORDER:
        cfg = CASES.get(name)
        if cfg is None or cfg["kind"] != "floor":
            continue
        files = sorted(glob.glob(os.path.join(DATA_DIR, cfg["pattern"])))
        files, dropped = unique_files(files)
        if not files:
            print(f"\n[{name}] no files matching {cfg['pattern']} - skipped")
            continue

        per_file = []
        spec = None
        for p in files:
            fl, sp = analyse_floor_file(p, cfg["f0"])
            per_file.append(fl)
            if spec is None:
                spec = sp
        # median across runs, so one bad capture cannot raise the floor
        merged = {}
        for h in range(2, N_HARMONICS + 1):
            vals = [d[h] for d in per_file if h in d]
            if vals:
                merged[h] = float(np.median(vals))
        floors[name] = merged
        floor_spectra[name] = spec

        print(f"\n[{name}] {cfg['note']}")
        print(f"  files: {len(files)}")
        print("  floor power at each harmonic frequency (arbitrary units):")
        for h, val in merged.items():
            print(f"    H{h} at {h*cfg['f0']/1000:5.1f} kHz : {val:.4e}")

    # ---- pass 2: signal cases -----------------------------------------------
    for name in RUN_ORDER:
        cfg = CASES.get(name)
        if cfg is None or cfg["kind"] != "signal":
            continue
        files = sorted(glob.glob(os.path.join(DATA_DIR, cfg["pattern"])))
        files, dropped = unique_files(files)
        if dropped:
            print(f"\n[{name}] dropped {len(dropped)} duplicate file(s): "
                  f"{', '.join(dropped)}")
        if not files:
            print(f"\n[{name}] no files matching {cfg['pattern']} - skipped")
            continue

        ext = floors.get(cfg["floor"]) if cfg["floor"] else None
        runs = []
        for p in files:
            try:
                runs.append(analyse_run(p, cfg["f0"], ext, cfg["lpf_corr"]))
            except Exception as exc:
                print(f"  !! {os.path.basename(p)}: {exc}")

        agg = aggregate(runs)
        if agg is None:
            continue
        results[name] = agg

        print("\n" + "-" * 74)
        print(f"[{name}]  f0 = {cfg['f0']:.0f} Hz, "
              f"{cfg['codes'] if cfg['codes'] else '-'} codes, "
              f"Point {cfg['point']}, {agg['n_runs']} runs")
        print(f"  {cfg['note']}")
        print(f"  scope rate {agg['fs']/1e3:.1f} kSa/s, {agg['n']} samples, "
              f"bin {agg['df']:.2f} Hz, cycles in record "
              f"{agg['f0']*agg['n']/agg['fs']:.1f}")
        print(f"  detected fundamental {agg['f0']:.2f} Hz, "
              f"Vpp {agg['vpp']:.3f} V")

        # ---- amplitude and clipping check ----------------------------------
        exp = expected_vpp(cfg["codes"])
        if exp is not None:
            err_pct = 100.0 * (agg["vpp"] - exp) / exp
            print(f"  expected Vpp {exp:.3f} V, measured {agg['vpp']:.3f} V "
                  f"({err_pct:+.2f} %)")
            if err_pct < -1.0:
                print("  *** WARNING: measured amplitude is more than 1 % "
                      "below theory.")
                print("      The signal may have run off the oscilloscope "
                      "screen. Recapture at a coarser volts/div.")

        cp, cn = agg["crest_pos"], agg["crest_neg"]
        print(f"  crest factor: positive {cp:.3f}, negative {cn:.3f} "
              f"(pure sine = 1.414)")
        if cp < 1.36 and cn < 1.36:
            print("  *** WARNING: both peaks are flattened - symmetric "
                  "limiting.")
            print("      This creates odd harmonics (H3, H5, H7) that are "
                  "not real distortion.")
        elif cp < 1.36:
            print("  *** WARNING: positive peaks flattened, negative peaks "
                  "are not.")
            print("      Asymmetric limiting, for example the MCP6002 input "
                  "handover near Vdd - 1.2 V.")
        elif cn < 1.36:
            print("  *** WARNING: negative peaks flattened, positive peaks "
                  "are not.")
            print("      Asymmetric limiting in the output stage.")

        if cfg["lpf_corr"]:
            print("  LPF response de-embedded from every component")

        print("  harmonic breakdown (mean over runs):")
        for h in range(2, N_HARMONICS + 1):
            r = agg["ratio"].get(h, 0.0)
            nm = agg["n_measured"].get(h, 0)
            tag = "measured" if nm > agg["n_runs"] / 2 else "at/below floor"
            print(f"    H{h} {h*agg['f0']/1000:7.2f} kHz  "
                  f"{dbc(r):8.2f} dBc   {tag} ({nm}/{agg['n_runs']} runs)")

        thd = agg["thd_f"]
        thd_r = thd / math.sqrt(1 + thd ** 2)
        print(f"  THD_F : {thd*100:.4f} %   ({db(thd):.2f} dB)")
        print(f"  THD_R : {thd_r*100:.4f} %")
        if agg["n_runs"] > 1:
            print(f"  spread over runs : +/- {agg['thd_db_std']:.2f} dB "
                  f"(1 sigma, n = {agg['n_runs']})")
        if agg["dominant"]:
            print(f"  dominant harmonic : H{agg['dominant']}")
            if agg["dominant"] == 2:
                print("    -> asymmetric distortion, output stage off centre "
                      "or near a rail")
            elif agg["dominant"] == 3:
                print("    -> symmetric compression, DAC INL curvature or "
                      "gentle limiting")

        n_lim = sum(1 for h in range(2, N_HARMONICS + 1)
                    if agg["n_measured"].get(h, 0) <= agg["n_runs"] / 2)
        if n_lim >= 3:
            print("  NOTE: most harmonics sit at or below the floor. "
                  "Report THD as an UPPER BOUND.")

    # ---- matched subset for the 5 kHz cases ---------------------------------
    for name in ("T3", "T4"):
        if name in results and results[name]["n_runs"] > MATCHED_SUBSET_N:
            sub = aggregate(results[name]["runs"][:MATCHED_SUBSET_N])
            print(f"\n[{name}] matched {MATCHED_SUBSET_N}-run subset "
                  f"(for fair comparison with the 1 kHz cases): "
                  f"THD = {sub['thd_f']*100:.4f} % ({db(sub['thd_f']):.2f} dB)")

    # ---- raw versus corrected at 5 kHz --------------------------------------
    if "T4" in results:
        cfg = dict(CASES["T4"])
        files = sorted(glob.glob(os.path.join(DATA_DIR, cfg["pattern"])))
        raw_runs = [analyse_run(p, cfg["f0"], None, False) for p in files]
        raw = aggregate(raw_runs)
        cor = results["T4"]
        print("\n" + "-" * 74)
        print("[T4] filter correction, Table 2")
        print(f"  THD from raw levels      : {raw['thd_f']*100:.4f} % "
              f"({db(raw['thd_f']):.2f} dB)")
        print(f"  THD after LPF correction : {cor['thd_f']*100:.4f} % "
              f"({db(cor['thd_f']):.2f} dB)")
        print(f"  difference               : "
              f"{db(cor['thd_f'])-db(raw['thd_f']):+.2f} dB")
        results["T4_raw"] = raw

    write_tables(results)
    make_figures(results, floors, floor_spectra)
    print("\nFigures written to ./" + FIG_DIR)
    print("Numbers written to results_thd.csv")


# =============================================================================
# 6.  OUTPUT TABLES
# =============================================================================

def write_tables(results):
    lines = ["case,point,f0_Hz,codes,n_runs,"
             + ",".join(f"H{h}_dBc" for h in range(2, N_HARMONICS + 1))
             + ",THD_percent,THD_dB,THD_std_dB,dominant"]

    for name in RUN_ORDER:
        if name not in results:
            continue
        cfg = CASES[name]
        a = results[name]
        hs = ",".join(f"{dbc(a['ratio'].get(h,0.0)):.2f}"
                      for h in range(2, N_HARMONICS + 1))
        lines.append(f"{name},{cfg['point']},{cfg['f0']:.0f},"
                     f"{cfg['codes'] if cfg['codes'] else ''},{a['n_runs']},"
                     f"{hs},{a['thd_f']*100:.4f},{db(a['thd_f']):.2f},"
                     f"{a['thd_db_std']:.2f},"
                     f"{'H'+str(a['dominant']) if a['dominant'] else ''}")

    with open("results_thd.csv", "w") as fh:
        fh.write("\n".join(lines) + "\n")

    # ---- Table 3, the amplitude sweep ---------------------------------------
    print("\n" + "=" * 74)
    print(" TABLE 3  --  amplitude sweep, f0 = 1 kHz")
    print("=" * 74)
    print(f"{'case':6} {'point':6} {'codes':6} {'swing (V)':14} "
          f"{'THD %':10} {'THD dB':9} {'dominant':9}")
    for grp in (SWEEP_B, SWEEP_A):
        for name in grp:
            if name not in results:
                continue
            cfg = CASES[name]
            a = results[name]
            c = cfg["codes"]
            lo = 2048 - c
            hi = 2048 + c
            swing = f"{lo*VDDA/DAC_FULLSCALE:.2f} - {hi*VDDA/DAC_FULLSCALE:.2f}"
            print(f"{name:6} {cfg['point']:6} {c:<6} {swing:14} "
                  f"{a['thd_f']*100:<10.4f} {db(a['thd_f']):<9.2f} "
                  f"{'H'+str(a['dominant']) if a['dominant'] else '-':9}")
        print()


# =============================================================================
# 7.  FIGURES
# =============================================================================

def _save(fig, stem):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG_DIR, f"{stem}.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close(fig)


def make_figures(results, floors, floor_spectra):

    # ---- Figure 1: spectrum of the headline case ----------------------------
    if "T2" in results:
        a = results["T2"]
        r = a["runs"][0]
        f, P = r["f_axis"], r["P_axis"]
        PdB = 10 * np.log10(P / np.max(P) + 1e-30)

        fig, ax = plt.subplots(figsize=(9, 5))
        ax.plot(f / 1e3, PdB, lw=0.7)
        for h in range(1, N_HARMONICS + 1):
            fx = h * a["f0"] / 1e3
            ax.axvline(fx, ls="--", alpha=0.35, lw=0.8)
            ax.text(fx, 4, f"H{h}", ha="center", fontsize=8)
        if "T9a" in floor_spectra and floor_spectra["T9a"]:
            sp = floor_spectra["T9a"]
            fl_dB = 10 * np.log10(sp["P_axis"] / np.max(P) + 1e-30)
            ax.plot(sp["f_axis"] / 1e3, fl_dB, lw=0.6, alpha=0.7,
                    label="noise floor (T9a)")
            ax.legend(fontsize=8)
        ax.set_xlabel("Frequency (kHz)")
        ax.set_ylabel("Magnitude (dB relative to fundamental)")
        ax.set_title("THD spectrum, Point B, 1 kHz, 1200 codes (T2)")
        ax.set_xlim(0, (N_HARMONICS + 1) * a["f0"] / 1e3)
        ax.set_ylim(-110, 10)
        ax.grid(alpha=0.3)
        _save(fig, "fig_thd_spectrum_T2")

    # ---- Figure 2: Point A versus Point B -----------------------------------
    if "T1" in results and "T2" in results:
        fig, ax = plt.subplots(figsize=(9, 5))
        for name, lab in (("T1", "Point A (raw DAC)"),
                          ("T2", "Point B (filter output)")):
            r = results[name]["runs"][0]
            PdB = 10 * np.log10(r["P_axis"] / np.max(r["P_axis"]) + 1e-30)
            ax.plot(r["f_axis"] / 1e3, PdB, lw=0.7, label=lab, alpha=0.85)
        ax.set_xlabel("Frequency (kHz)")
        ax.set_ylabel("Magnitude (dB relative to fundamental)")
        ax.set_title("Spectrum before and after the low-pass filter, 1 kHz")
        ax.set_xlim(0, 10)
        ax.set_ylim(-110, 10)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=9)
        _save(fig, "fig_thd_A_vs_B")

    # ---- Figure 3: harmonic bar chart ---------------------------------------
    if "T2" in results:
        a = results["T2"]
        hs = list(range(2, N_HARMONICS + 1))
        vals = [dbc(a["ratio"].get(h, 0.0)) for h in hs]
        meas = [a["n_measured"].get(h, 0) > a["n_runs"] / 2 for h in hs]
        fig, ax = plt.subplots(figsize=(7, 4.5))
        bars = ax.bar([f"H{h}" for h in hs], vals)
        for b, m in zip(bars, meas):
            if not m:
                b.set_hatch("//")
                b.set_alpha(0.55)
        ax.set_ylabel("Level (dBc)")
        ax.set_title("Individual harmonic levels, Point B, 1 kHz, 1200 codes\n"
                     "(hatched = at or below the noise floor, upper bound)")
        ax.grid(alpha=0.3, axis="y")
        _save(fig, "fig_thd_harmonics_T2")

    # ---- Figure 4: THD versus amplitude -------------------------------------
    have_b = [n for n in SWEEP_B if n in results]
    have_a = [n for n in SWEEP_A if n in results]
    if have_b:
        fig, ax = plt.subplots(figsize=(8, 5))
        for grp, lab, mk in ((have_b, "Point B (filter output)", "o"),
                             (have_a, "Point A (raw DAC)", "s")):
            if not grp:
                continue
            xs = [CASES[n]["codes"] for n in grp]
            ys = [db(results[n]["thd_f"]) for n in grp]
            es = [results[n]["thd_db_std"] for n in grp]
            order = np.argsort(xs)
            xs = np.array(xs)[order]
            ys = np.array(ys)[order]
            es = np.array(es)[order]
            ax.errorbar(xs, ys, yerr=es, marker=mk, capsize=3, label=lab)
        ax.axvline(1200, ls=":", alpha=0.6)
        ax.text(1210, ax.get_ylim()[0] + 1, "PM operating point", fontsize=8)
        ax.set_xlabel("DAC amplitude (codes)")
        ax.set_ylabel("THD (dB relative to fundamental)")
        ax.set_title("THD versus DAC amplitude, f0 = 1 kHz\n"
                     "(more negative is better)")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=9)
        _save(fig, "fig_thd_vs_amplitude")

    # ---- Figure 5: 5 kHz raw versus corrected -------------------------------
    if "T4" in results and "T4_raw" in results:
        hs = list(range(2, N_HARMONICS + 1))
        raw = [dbc(results["T4_raw"]["ratio"].get(h, 0.0)) for h in hs]
        cor = [dbc(results["T4"]["ratio"].get(h, 0.0)) for h in hs]
        x = np.arange(len(hs))
        fig, ax = plt.subplots(figsize=(7.5, 4.5))
        ax.bar(x - 0.2, raw, 0.4, label="raw")
        ax.bar(x + 0.2, cor, 0.4, label="filter corrected")
        ax.set_xticks(x)
        ax.set_xticklabels([f"H{h}" for h in hs])
        ax.set_ylabel("Level (dBc)")
        ax.set_title("Harmonic levels at 5 kHz, before and after "
                     "de-embedding the filter (T4)")
        ax.grid(alpha=0.3, axis="y")
        ax.legend(fontsize=9)
        _save(fig, "fig_thd_5kHz_correction")


if __name__ == "__main__":
    main()
