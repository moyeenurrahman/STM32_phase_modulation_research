"""
Frequency accuracy analysis - reads Tektronix TBS1000C CSV exports.

Usage:
    python analyse_freq.py 5000 TEK0000CH1.CSV [TEK0000CH2.CSV]

The first argument is the target frequency in Hz. Giving two files runs
the T7 comparison between the filtered output and the raw DAC output.
"""

import os
import sys
import numpy as np

SCOPE_TIMEBASE_PPM = 20.0     # Tektronix TBS1000C specification
RAW_DAC_TOLERANCE_PPM = 60.0  # demonstrated spread when the zero-crossing
                              # estimator is applied to the raw DAC
                              # staircase at ~25 samples per carrier cycle


# ============================================================
# Loader - tolerant of the several layouts Tektronix scopes use
# ============================================================
def load_tek_csv(path):
    """Read a scope CSV. Copes with several Tektronix / OpenChoice layouts:
    header block in columns 1-2, plain two-column time+volts, volts-only
    single column, and comma, semicolon or tab delimiters."""
    with open(path, "r", errors="ignore") as f:
        lines = [ln.rstrip("\n").rstrip("\r") for ln in f]
    lines = [ln for ln in lines if ln.strip() != ""]
    if not lines:
        raise ValueError(f"{path}: file is empty")

    # --- pick the delimiter by which one splits the file most consistently
    sample = lines[max(0, len(lines) // 2): len(lines) // 2 + 40]
    best, delim = -1, ","
    for d in (",", ";", "\t", " "):
        counts = [ln.count(d) for ln in sample]
        if counts and max(counts) > 0:
            score = sum(1 for c in counts if c == max(counts))
            if score > best:
                best, delim = score, d

    rows = [[p.strip().strip('"') for p in ln.split(delim)] for ln in lines]

    header = {}
    for r in rows:
        if r and r[0] and not _is_float(r[0]) and len(r) > 1 and r[1]:
            header.setdefault(r[0], r[1])

    width = max(len(r) for r in rows)
    grid = np.full((len(rows), width), np.nan)
    for i, r in enumerate(rows):
        for j, p in enumerate(r):
            p = p.replace(" ", "")
            if _is_float(p):
                grid[i, j] = float(p)

    # exact sample interval from the header, if the scope wrote one
    dt = None
    for key in header:
        if key.lower().replace(" ", "") in (
                "sampleinterval", "sampleperiod", "horizontalinterval",
                "xincr", "x_increment", "xincrement"):
            if _is_float(header[key]):
                dt = float(header[key])
                break

    # --- find columns that are numeric in most of the file
    frac = np.mean(~np.isnan(grid), axis=0)
    cand = [c for c in range(width) if frac[c] > 0.8]
    cand = [c for c in cand
            if np.nanstd(grid[:, c]) > 0 and not np.all(np.isnan(grid[:, c]))]

    if not cand:
        _dump_head(path, lines)
        raise ValueError(f"{path}: no numeric data columns found")

    # Drop rows that do not have every data column filled. A header field
    # containing a comma shifts one row sideways, and a single shifted row
    # is enough to make the time column look non-uniform.
    rowmask = np.all(~np.isnan(grid[:, cand]), axis=1)
    grid = grid[rowmask]

    # a time column rises steadily with an almost constant step
    t_col = None
    for c in cand:
        col = grid[:, c]
        if len(col) < 10:
            continue
        d = np.diff(col)
        step = np.median(d)
        # robust: most steps close to the median, rather than every step
        good = np.mean(np.abs(d - step) < 0.2 * abs(step))
        if (step > 0 and good > 0.95 and np.all(d > 0)
                and 1e-12 < step < 1.0):
            t_col = c
            break

    if t_col is None:
        # volts-only file: rebuild the time axis from the sample interval
        v_col = max(cand, key=lambda c: np.nanstd(grid[:, c]))
        v = grid[~np.isnan(grid[:, v_col]), v_col]
        if dt is None:
            _dump_head(path, lines)
            raise ValueError(
                f"{path}: only one data column and no sample interval in the "
                f"header. Set MANUAL_DT at the bottom of this script.")
        t = dt * np.arange(len(v))
        return dict(t=t, v=v, channels=[v], dt=dt, header=header, path=path,
                    delim=delim, note="time axis rebuilt from sample interval")

    others = [c for c in cand if c != t_col]
    if not others:
        _dump_head(path, lines)
        raise ValueError(f"{path}: found a time column but no voltage column")

    # a real waveform column swings; filler columns are flat
    span = {c: np.nanmax(grid[:, c]) - np.nanmin(grid[:, c]) for c in others}
    biggest = max(span.values())
    chans = sorted(c for c in others if span[c] > 0.05 * biggest)

    mask = ~np.isnan(grid[:, t_col])
    for c in chans:
        mask &= ~np.isnan(grid[:, c])

    vs = [grid[mask, c] for c in chans]
    return dict(t=grid[mask, t_col], v=vs[0], channels=vs, dt=dt,
                header=header, path=path, delim=delim, note="")


def _dump_head(path, lines):
    print(f"\n--- first 12 lines of {path} ---")
    for ln in lines[:12]:
        print("   " + ln[:120])
    print("--- end ---\n")


def _is_float(s):
    try:
        float(s)
        return True
    except (TypeError, ValueError):
        return False


# ============================================================
# Hardened zero-crossing estimator (validated in Session 0)
# ============================================================
def is_sample_and_hold(v, dt, f_set):
    """Detect a raw DAC staircase.

    A reconstructed sine changes by many LSBs between scope samples near
    the zero crossing. A sample-and-hold waveform instead sits on one code
    for a whole DAC period, so consecutive scope samples repeat. When the
    scope rate and the DAC rate are commensurate those flat plateaus snap
    the interpolated crossings onto a fixed grid, which locks the fitted
    frequency near its nominal value and makes a badly biased result look
    clean. The residual-based sigma is meaningless in that case.
    """
    v = np.asarray(v, float)
    mid = np.abs(v - np.mean(v)) < 0.3 * np.std(v)   # near the crossings
    if mid.sum() < 20:
        return False
    pairs = mid[:-1] & mid[1:]          # both samples near the crossing
    if pairs.sum() < 20:
        return False
    same = np.mean(np.abs(np.diff(v)[pairs]) < 1e-9)
    return same > 0.15


def freq_zero_crossing(t, v, dt=None, hyst=0.25, trim=0.01):
    """trim discards this fraction of samples from each end of the record.
    Scopes commonly corrupt the first and last few samples, and a single
    bad crossing there shifts the cycle numbering for the whole fit."""
    v = np.asarray(v, float)
    t = np.asarray(t, float)
    if dt is not None:
        t = t[0] + dt * np.arange(len(t))
    # Scopes can freeze on a constant value at the start or end of a
    # record. Find those runs and cut them, then take the fixed trim.
    def stuck_run(x):
        if len(x) < 10:
            return 0
        tol = 1e-9 + 0.002 * (np.max(x) - np.min(x))
        n = 0
        while n + 1 < len(x) and abs(x[n + 1] - x[0]) <= tol:
            n += 1
        return n + 1 if n >= 4 else 0

    lead = stuck_run(v)
    tail = stuck_run(v[::-1])
    if lead or tail:
        v = v[lead:len(v) - tail]
        t = t[lead:len(t) - tail]

    if trim > 0:
        n = max(3, int(round(trim * len(v))))
        v = v[n:-n]
        t = t[n:-n]
    v = v - np.mean(v)

    thr = hyst * np.sqrt(2.0) * np.std(v)
    raw = np.where((v[:-1] < 0) & (v[1:] >= 0))[0]
    if len(raw) < 10:
        raise ValueError("Too few zero crossings")

    keep, last = [], -1
    for i in raw:
        # window starts AFTER the previous accepted crossing, otherwise the
        # negative sample belonging to that crossing re-arms the trigger and
        # a glitch dipping back through zero is wrongly accepted
        seg = v[last + 1:i + 1]
        if seg.size and seg.min() < -thr:
            keep.append(i)
            last = i
    keep = np.array(keep)

    t0 = t[keep] - v[keep] * (t[keep + 1] - t[keep]) / (v[keep + 1] - v[keep])

    # Cycle index from time, not from counting. The starting period estimate
    # is the record average rather than the median: on a staircase waveform
    # the individual periods are skewed and the median is biased. Two
    # refinement passes remove any residual bias.
    T_nom = (t0[-1] - t0[0]) / (len(t0) - 1)
    for _ in range(3):
        k = np.round((t0 - t0[0]) / T_nom)
        A = np.vstack([k, np.ones_like(k)]).T
        coef, *_ = np.linalg.lstsq(A, t0, rcond=None)
        T_nom = coef[0]

    period = coef[0]
    res = t0 - (A @ coef)

    ok = True
    notes = []
    if np.max(np.abs(res)) > 0.25 * period:
        ok, _ = False, notes.append("residual > quarter period")
    if len(np.unique(k)) != len(k):
        ok, _ = False, notes.append("duplicate cycle indices")

    # Precision of the fitted slope, from the standard error of a
    # least-squares line through M equally spaced points.
    rms = np.sqrt(np.mean(res ** 2))
    M = len(t0)
    sigma_ppm = rms * np.sqrt(12.0) / (M ** 1.5) / period * 1e6

    return dict(freq=1.0 / period, jitter_rms=rms, sigma_ppm=sigma_ppm,
                n_cycles=M, span=t0[-1] - t0[0], ok=ok,
                notes="; ".join(notes) if notes else "clean")


def freq_fft_interp(t, v, dt=None):
    """METHOD C - interpolated FFT peak.

    Used when modulation is active. Phase modulation displaces the zero
    crossings, which biases the regression, but it leaves the carrier
    spectral line where it is. A Blackman window keeps the side lobes low
    enough that the message sidebands do not pull the peak, and parabolic
    interpolation on the log magnitudes locates it to a small fraction
    of a bin.
    """
    v = np.asarray(v, float) - np.mean(v)
    if dt is None:
        dt = (t[-1] - t[0]) / (len(t) - 1)
    N = len(v)
    V = np.abs(np.fft.rfft(v * np.blackman(N)))
    V[0] = 0.0
    k = int(np.argmax(V))
    if k < 1 or k >= len(V) - 1:
        raise ValueError("carrier peak sits at the edge of the spectrum")
    a, b, c = (np.log(max(V[k - 1], 1e-30)), np.log(max(V[k], 1e-30)),
               np.log(max(V[k + 1], 1e-30)))
    delta = 0.5 * (a - c) / (a - 2 * b + c)
    freq = (k + delta) / (N * dt)

    # sideband check, so a modulated capture is never mistaken for a clean one
    order = np.argsort(V)[::-1]
    second = next((i for i in order if abs(i - k) > 2), k)
    sb_dB = 20 * np.log10(V[second] / V[k])

    return dict(freq=freq, bin_hz=1.0 / (N * dt), n_cycles=freq * N * dt,
                jitter_rms=0.0, sigma_ppm=float("nan"),
                strongest_other_dB=sb_dB, ok=True, span=N * dt,
                notes=f"strongest other component {sb_dB:.1f} dB")


def ppm(f_meas, f_set):
    return (f_meas - f_set) / f_set * 1e6


# ============================================================
def analyse(path, f_set, ch=0, d=None):
    if d is None:
        d = load_tek_csv(path)
    if d["dt"] is None and MANUAL_DT:
        d["dt"] = MANUAL_DT
    d = dict(d, v=d["channels"][ch])
    src = d["header"].get("Source", "?")
    label = f"channel {ch + 1}" if len(d["channels"]) > 1 else f"Source: {src}"
    print(f"\n--- {d['path']}   ({label}) ---")
    print(f"  points          : {len(d['t'])}")
    if d["dt"] is None:
        print("  sample interval : NOT in header, using time column")
    else:
        print(f"  sample interval : {d['dt']*1e9:.1f} ns  (from header)")
        print(f"  scope rate      : {1/d['dt']/1e3:.2f} kSa/s")
    print(f"  Vpp             : {d['v'].max()-d['v'].min():.3f} V")

    span = d["t"][-1] - d["t"][0]
    problems = []
    if len(d["t"]) < 500:
        problems.append(f"only {len(d['t'])} points - a scope record should "
                        f"have thousands")
    if f_set * span < 20:
        problems.append(f"record holds only {f_set*span:.1f} cycles of the "
                        f"target frequency - use a slower timebase")
    if not (1e-7 < span < 100.0):
        problems.append(f"record span is {span:.3e} s, which is not a "
                        f"plausible capture length")
    if problems:
        print("  FILE DOES NOT LOOK LIKE A SCOPE CAPTURE:")
        for m in problems:
            print("    - " + m)
        print("  Nothing below this line is meaningful.")

    if PM_MODE:
        r = freq_fft_interp(d["t"], d["v"], dt=d["dt"])
        print(f"  method          : C, interpolated FFT (modulation active)")
        print(f"  FFT bin spacing : {r['bin_hz']:.2f} Hz")
        print(f"  record span     : {r['span']*1000:.2f} ms, "
              f"{r['n_cycles']:.0f} cycles")
        print(f"  carrier         : {r['freq']:.4f} Hz")
        print(f"  ERROR           : {ppm(r['freq'], f_set):+.2f} ppm")
        print(f"  strongest sideband: {r['strongest_other_dB']:.1f} dB")
        if r["strongest_other_dB"] < -20:
            print("  >> No strong sidebands. Is KP really 1.0f?")
    else:
        r = freq_zero_crossing(d["t"], d["v"], dt=d["dt"])
        print(f"  method          : D, zero-crossing regression")
        print(f"  record span     : {r['span']*1000:.2f} ms, "
              f"{r['n_cycles']} cycles")
        print(f"  frequency       : {r['freq']:.4f} Hz")
        print(f"  ERROR           : {ppm(r['freq'], f_set):+.2f} ppm")
        print(f"  crossing jitter : {r['jitter_rms']*1e9:.0f} ns")
        r["sah"] = is_sample_and_hold(d["v"], d["dt"], f_set)
        if r["sah"]:
            print(f"  estimator 1 sigma: {r['sigma_ppm']:.3f} ppm "
                  f"(NOT VALID - see below)")
            print("  >> SAMPLE-AND-HOLD WAVEFORM (raw DAC).")
            print("     The staircase sits on one code for a full DAC period")
            print("     at each zero crossing, so the crossings snap to the")
            print("     scope sample grid. The quoted sigma is far too")
            print("     optimistic; the true uncertainty is tens of ppm.")
            print("     Use the filtered channel for the reported result.")
        else:
            print(f"  estimator 1 sigma: {r['sigma_ppm']:.3f} ppm")
    print(f"  status          : {'OK' if r['ok'] else 'CHECK'}  ({r['notes']})")

    err = abs(ppm(r["freq"], f_set))
    if not PM_MODE and r["jitter_rms"] > 0.02 / f_set:
        pass
    if err > 100000 or problems:
        print("  >> RESULT REJECTED. The measured frequency is nowhere near")
        print("     the target, so the file or the target setting is wrong.")
        r["ok"] = False
    elif not r["ok"]:
        print("  >> Do not use this result. Recapture.")
    if not PM_MODE and r["jitter_rms"] > 0.02 / f_set:
        print("  >> Jitter is high. Check probe, coupling, vertical scale.")
    return r


# ============================================================
# LAB MODE - edit these three lines, then press Run in PyCharm.
# Command-line arguments, if given, override them.
# ============================================================
F_SET_DEFAULT = 5000.0
FILE_1 = "capture_CH1.csv"     # filtered output
FILE_2 = ""                    # raw DAC, leave "" if only one channel
MANUAL_DT = None               # e.g. 7.5e-6, only if the CSV has no header
                               # and no time column
PM_MODE = False                # True for T6 only. Switches to Method C,
                               # because modulation displaces the zero
                               # crossings and biases Method D by tens of ppm


if __name__ == "__main__":
    if len(sys.argv) >= 3:
        f_set = float(sys.argv[1])
        files = sys.argv[2:]
    else:
        f_set = F_SET_DEFAULT
        files = [FILE_1] + ([FILE_2] if FILE_2 else [])
        print("(lab mode - using the settings at the bottom of the file)")

    missing = [f for f in files if not os.path.exists(f)]
    if missing:
        print("\nFile not found: " + ", ".join(missing))
        print("Working folder: " + os.getcwd())
        print("CSV files here: " +
              (", ".join(f for f in os.listdir(".") if f.lower().endswith(".csv"))
               or "none"))
        sys.exit(1)

    print(f"Target frequency: {f_set:.3f} Hz")

    if len(files) == 1:
        d = load_tek_csv(files[0])
        n = len(d["channels"])
        if n > 1:
            print(f"({n} channels found in this file)")
        results = [analyse(files[0], f_set, ch=i, d=d) for i in range(n)]
    else:
        results = [analyse(p, f_set) for p in files]

    if len(results) == 2:
        a, b = results
        gap = ppm(a["freq"], f_set) - ppm(b["freq"], f_set)
        print("\n=== T7 verdict ===")
        print(f"  CH1 (filtered) : {ppm(a['freq'], f_set):+.2f} ppm  "
              f"(1 sigma {a['sigma_ppm']:.3f})")
        print(f"  CH2 (raw DAC)  : {ppm(b['freq'], f_set):+.2f} ppm  "
              f"(1 sigma {b['sigma_ppm']:.3f})")
        print(f"  difference     : {gap:+.2f} ppm")

        # If either channel is a staircase, the comparison is limited by how
        # well that staircase can be measured, not by the scope timebase.
        sah = any(x.get("sah") for x in results)
        tol = RAW_DAC_TOLERANCE_PPM if sah else SCOPE_TIMEBASE_PPM
        print(f"  tolerance      : +/- {tol:.0f} ppm "
              + ("(set by the staircase channel)" if sah
                 else "(scope timebase)"))
        if abs(gap) < tol:
            print("  PASS - the two channels agree within the uncertainty of")
            print("         the noisier one, so no evidence that the filter")
            print("         shifts the frequency.")
        else:
            print("  FAIL - investigate before continuing.")
        print("\n  Note: the raw DAC channel is expected to show larger")
        print("  crossing jitter. At 20 samples per carrier cycle the DAC")
        print("  sits exactly on mid-code (2048) for one full sample period")
        print("  at each zero crossing, so the crossing instant is uncertain")
        print("  by up to one 10 us hold. This is why the filtered channel is")
        print("  the one reported.")
