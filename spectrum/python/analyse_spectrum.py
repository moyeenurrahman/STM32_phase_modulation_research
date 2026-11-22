#!/usr/bin/env python3
"""
analyse_spectrum.py
-------------------
Spectral analysis for the baseband PM thesis.

Handles both capture families:

    Family W : 400 us/div -> 312.5 kSa/s, 156.25 Hz bins, Nyquist 156.25 kHz
               Job: DAC zero-order-hold images at 95 and 105 kHz,
                    before/after LPF comparison.

    Family N : 2 ms/div   ->  62.5 kSa/s,  31.25 Hz bins, Nyquist  31.25 kHz
               Job: PM sidebands, Bessel carrier null, B99%, harmonics,
                    BPSK main lobe and sidelobes.

Records are 2000 points, so a single periodogram over the full record is
used. Welch segmentation is impossible at this record length without
destroying the frequency resolution.

WARNING - Family N aliasing:
    At 62.5 kSa/s the 105 kHz image folds to exactly 20 kHz and the 95 kHz
    image folds to 30 kHz. At Point B these are suppressed and harmless.
    NEVER analyse a Point A capture in Family N - the folded image would
    sit on top of the 4th harmonic.

System parameters assumed:
    fc = 5 kHz, fm = 500 Hz, Rs = 500 bps
    DAC fs = 100 kSa/s, buffer N = 1000  ->  pattern repeats at 100 Hz
"""

import hashlib
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import welch
from scipy.special import jv


# ======================================================================
# CONFIGURATION
# ======================================================================

DATA_DIR = "."

# Define every case here. 'pattern' is a substring matched against
# filenames. 'kind' selects which analysis runs.
#
#   kind = "tone"  unmodulated carrier
#   kind = "pm"    analog phase modulation
#   kind = "bpsk"  phase-shift keying
#   kind = "floor" output disabled, instrument noise floor
#
CASES = {
    # ---- Family W (already captured during the SNR session) ----
    "S1":  dict(pattern="S1",  family="W", kind="tone",  point="A", kp=0.0),
    "S2":  dict(pattern="S2",  family="W", kind="tone",  point="B", kp=0.0),
    "S6W": dict(pattern="N6",  family="W", kind="bpsk",  point="B", kp=None),
    "S8W": dict(pattern="S8W", family="W", kind="floor", point="B", kp=None),

    # ---- Family N (to be captured) ----
    "S2N": dict(pattern="S2N", family="N", kind="tone",  point="B", kp=0.0),
    "S3":  dict(pattern="S3",  family="N", kind="pm",    point="B", kp=1.0),
    "S4":  dict(pattern="S4",  family="N", kind="pm",    point="B", kp=2.4048),
    "S5":  dict(pattern="S5",  family="N", kind="pm",    point="B", kp=3.8317),
    "S6N": dict(pattern="S6N", family="N", kind="bpsk",  point="B", kp=None),
    "S7N": dict(pattern="S7N", family="N", kind="bpsk",  point="B", kp=None),
    "S8N": dict(pattern="S8N", family="N", kind="floor", point="B", kp=None),
}

# Which floor case belongs to which family
FLOOR_FOR_FAMILY = {"W": "S8W", "N": "S8N"}

FAMILY_FS = {"W": 312_500.0, "N": 62_500.0}
FS_TOLERANCE = 0.01                 # 1 percent

FC = 5_000.0                        # carrier
FM = 500.0                          # modulating frequency
RS = 500.0                          # symbol rate

N_SIDEBANDS = 6                     # analyse fc +/- n*fm for n = 1..6
N_HARMONICS = 6                     # 2fc .. 6fc

PEAK_TOL_BINS = 3                   # main-lobe half-width for peak power

B99_LOW = 2_000.0                   # symmetric about fc, excludes 2nd harmonic
B99_HIGH = 8_000.0

SHOULDER_OFFSET = 2_000.0           # measure at fc +/- this

IMAGE_LOW = 95_000.0
IMAGE_HIGH = 105_000.0

FLOOR_MARGIN_DB = 3.0               # below this, results are instrument-limited

SAVE_PLOTS = True
PLOT_DPI = 300


# ======================================================================
# FILE HANDLING
# ======================================================================

def md5_of(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def load_scope_csv(path):
    """
    Read a Tektronix TBS1000C CSV.

    Copes with the pendrive export (metadata block on the left, data
    columns further right) and the OpenChoice two-column export.
    Returns (t, v).
    """
    rows = []
    with open(path, "r", errors="ignore") as fh:
        for line in fh:
            out = []
            for p in line.split(","):
                try:
                    out.append(float(p.strip()))
                except ValueError:
                    out.append(np.nan)
            rows.append(out)

    if not rows:
        raise ValueError("file is empty")

    width = max(len(r) for r in rows)
    arr = np.full((len(rows), width), np.nan)
    for i, r in enumerate(rows):
        arr[i, :len(r)] = r

    good = [c for c in range(width) if np.mean(np.isfinite(arr[:, c])) > 0.9]
    if len(good) < 2:
        raise ValueError("could not find two numeric columns")

    # Time column increases monotonically.
    t_col = None
    for c in good:
        col = arr[np.isfinite(arr[:, c]), c]
        if len(col) > 10 and np.all(np.diff(col) > 0):
            t_col = c
            break
    if t_col is None:
        t_col = good[-2]

    v_col = [c for c in good if c != t_col][-1]

    m = np.isfinite(arr[:, t_col]) & np.isfinite(arr[:, v_col])
    t, v = arr[m, t_col], arr[m, v_col]

    if len(t) < 100:
        raise ValueError(f"only {len(t)} usable rows")

    return t, v


def collect_case(case_id, cfg):
    """
    Load, screen and average every record belonging to one case.
    Returns dict with f, Pxx, n_used, n_dup, n_files, or None.
    """
    names = sorted(f for f in os.listdir(DATA_DIR)
                   if cfg["pattern"].lower() in f.lower()
                   and f.lower().endswith(".csv"))
    if not names:
        return None

    fs_want = FAMILY_FS[cfg["family"]]

    seen, kept, n_dup = {}, [], 0
    for n in names:
        p = os.path.join(DATA_DIR, n)
        h = md5_of(p)
        if h in seen:
            n_dup += 1
        else:
            seen[h] = p
            kept.append(p)

    spectra, n_bad = [], 0
    f_ref = None
    for p in kept:
        try:
            t, v = load_scope_csv(p)
        except Exception:
            n_bad += 1
            continue

        fs = 1.0 / np.mean(np.diff(t))
        if abs(fs - fs_want) / fs_want > FS_TOLERANCE:
            n_bad += 1
            continue
        if len(v) < 1500:
            n_bad += 1
            continue

        x = v - np.mean(v)
        f, Pxx = welch(x, fs=fs, window="blackmanharris",
                       nperseg=len(x), noverlap=0,
                       scaling="density", detrend=False)

        if f_ref is None:
            f_ref = f
        elif len(f) != len(f_ref):
            n_bad += 1
            continue

        spectra.append(Pxx)

    if not spectra:
        return None

    return dict(
        f=f_ref,
        Pxx=np.mean(spectra, axis=0),      # average power, not dB
        n_files=len(names),
        n_dup=n_dup,
        n_bad=n_bad,
        n_used=len(spectra),
        cfg=cfg,
    )


# ======================================================================
# SPECTRAL PRIMITIVES
# ======================================================================

def to_db(p):
    return 10.0 * np.log10(np.maximum(np.asarray(p, dtype=float), 1e-30))


def lobe_power(f, Pxx, f0, tol_bins=PEAK_TOL_BINS):
    """
    Power in the main lobe around f0, summed across +/- tol_bins.
    Uses summation rather than np.max, matching the THD analysis.
    """
    df = f[1] - f[0]
    m = np.abs(f - f0) <= tol_bins * df
    if not np.any(m):
        return 0.0
    return float(np.sum(Pxx[m]) * df)


def lobe_amplitude(f, Pxx, f0, tol_bins=PEAK_TOL_BINS):
    """Sinusoid amplitude from lobe power: P = A^2 / 2."""
    return np.sqrt(2.0 * lobe_power(f, Pxx, f0, tol_bins))


def band_power(f, Pxx, lo, hi):
    df = f[1] - f[0]
    m = (f >= lo) & (f <= hi)
    if not np.any(m):
        raise ValueError(f"no bins in {lo}-{hi} Hz")
    return float(np.sum(Pxx[m]) * df), f[m][0], f[m][-1], int(m.sum())


def occupied_bandwidth(f, Pxx, lo, hi, frac=0.99):
    """B99% inside a restricted band. Returns (B, f_low, f_high)."""
    df = f[1] - f[0]
    m = (f >= lo) & (f <= hi)
    fb, Pb = f[m], Pxx[m]
    cum = np.cumsum(Pb * df)
    total = cum[-1]
    tail = (1.0 - frac) / 2.0
    i_lo = int(np.searchsorted(cum, tail * total))
    i_hi = int(np.searchsorted(cum, (1.0 - tail) * total))
    i_hi = min(i_hi, len(fb) - 1)
    return fb[i_hi] - fb[i_lo], fb[i_lo], fb[i_hi]


def recover_kp(amps):
    """
    Recover the phase deviation index from the measured sideband pattern.

    amps[n] is the amplitude at fc +/- n*fm, n = 0..N.
    Theory: amplitude of order n is proportional to |J_n(kp)|.

    A least-squares fit over all orders is used rather than the simple
    J1/J0 ratio, because that ratio blows up at kp = 2.4048 where the
    carrier nulls. Fitting the whole pattern stays well behaved there.
    """
    a = np.asarray(amps, dtype=float)
    if np.all(a <= 0):
        return np.nan, np.nan
    a = a / np.linalg.norm(a)

    grid = np.arange(0.01, 6.0, 0.001)
    best_kp, best_err = np.nan, np.inf
    orders = np.arange(len(a))
    for kp in grid:
        th = np.abs(jv(orders, kp))
        nrm = np.linalg.norm(th)
        if nrm <= 0:
            continue
        err = float(np.sum((a - th / nrm) ** 2))
        if err < best_err:
            best_err, best_kp = err, kp
    return best_kp, np.sqrt(best_err)


def find_null(f, Pxx, f_target, search_hz=250.0):
    """Locate a local minimum near f_target. Used for BPSK main-lobe nulls."""
    m = np.abs(f - f_target) <= search_hz
    if not np.any(m):
        return np.nan, np.nan
    fb, Pb = f[m], Pxx[m]
    i = int(np.argmin(Pb))
    return float(fb[i]), float(Pb[i])


def find_peak(f, Pxx, lo, hi):
    m = (f >= lo) & (f <= hi)
    if not np.any(m):
        return np.nan, np.nan
    fb, Pb = f[m], Pxx[m]
    i = int(np.argmax(Pb))
    return float(fb[i]), float(Pb[i])


# ======================================================================
# PER-CASE ANALYSIS
# ======================================================================

def analyse_family_w(name, res, floor):
    f, P = res["f"], res["Pxx"]
    cfg = res["cfg"]

    print(f"\n  -- DAC images (Family W) --")

    if cfg["kind"] == "bpsk":
        p_sig, _, _, _ = band_power(f, P, FC - RS, FC + RS)
        ref_label = "main lobe"
    else:
        p_sig = lobe_power(f, P, FC)
        ref_label = "carrier"

    print(f"  {ref_label:<12}          : {to_db(p_sig):8.2f} dB(V^2)")

    for label, ftgt in (("Image  95 kHz", IMAGE_LOW),
                        ("Image 105 kHz", IMAGE_HIGH)):
        p_img = lobe_power(f, P, ftgt)
        dbc = to_db(p_img) - to_db(p_sig)

        note = ""
        if floor is not None:
            p_flr = lobe_power(floor["f"], floor["Pxx"], ftgt)
            margin = to_db(p_img) - to_db(p_flr)
            if margin < FLOOR_MARGIN_DB:
                note = f"  <- within {margin:.1f} dB of floor: LOWER BOUND"
            else:
                note = f"  ({margin:+.1f} dB above floor)"

        print(f"  {label}          : {to_db(p_img):8.2f} dB(V^2)  "
              f"{dbc:+8.2f} dBc{note}")

    if floor is None:
        print("  ** No noise floor loaded - suppression cannot be validated **")


def analyse_family_n(name, res, floor):
    f, P = res["f"], res["Pxx"]
    cfg = res["cfg"]
    kind = cfg["kind"]

    if cfg["point"] == "A":
        print("\n  ** REFUSING Family N analysis of a Point A capture. **")
        print("     The 105 kHz image folds to 20 kHz and corrupts H4.")
        return

    # ---------- harmonics (all kinds) ----------
    p_carrier = lobe_power(f, P, FC)
    print(f"\n  -- Harmonics --")
    for n in range(2, N_HARMONICS + 1):
        fh = n * FC
        if fh > f[-1]:
            break
        ph = lobe_power(f, P, fh)
        print(f"  H{n} at {fh/1e3:5.1f} kHz     : "
              f"{to_db(ph):8.2f} dB(V^2)  {to_db(ph)-to_db(p_carrier):+8.2f} dBc")

    # ---------- PM / tone ----------
    if kind in ("pm", "tone"):
        print(f"\n  -- Sidebands --")
        amps = [lobe_amplitude(f, P, FC)]
        print(f"  n=0  carrier      : {to_db(lobe_power(f,P,FC)):8.2f} dB(V^2)")
        for n in range(1, N_SIDEBANDS + 1):
            a_lo = lobe_amplitude(f, P, FC - n * FM)
            a_hi = lobe_amplitude(f, P, FC + n * FM)
            a = 0.5 * (a_lo + a_hi)
            amps.append(a)
            p_lo = lobe_power(f, P, FC - n * FM)
            p_hi = lobe_power(f, P, FC + n * FM)
            print(f"  n={n}  {FC-n*FM:5.0f}/{FC+n*FM:5.0f} Hz : "
                  f"{to_db(p_lo):8.2f} / {to_db(p_hi):8.2f} dB(V^2)   "
                  f"{to_db(p_lo)-to_db(p_carrier):+7.2f} / "
                  f"{to_db(p_hi)-to_db(p_carrier):+7.2f} dBc")

        if kind == "pm":
            kp_fit, resid = recover_kp(amps)
            kp_cmd = cfg["kp"]
            print(f"\n  -- Bessel recovery --")
            print(f"  Commanded kp      : {kp_cmd:.4f} rad")
            print(f"  Recovered kp      : {kp_fit:.4f} rad   "
                  f"(fit residual {resid:.4f})")
            if np.isfinite(kp_fit) and kp_cmd > 0:
                err = 100.0 * (kp_fit - kp_cmd) / kp_cmd
                flag = "" if abs(err) <= 10 else "   <-- OUTSIDE 10% TOLERANCE"
                print(f"  Error             : {err:+.2f} %{flag}")

            if abs(kp_cmd - 2.4048) < 0.01:
                depth = to_db(p_carrier) - to_db(np.median(P))
                print(f"  Carrier above median floor: {depth:.2f} dB")
                print("  (A deep null here is the headline result.)")

        # ---------- B99 ----------
        try:
            b99, fl, fh = occupied_bandwidth(f, P, B99_LOW, B99_HIGH)
            print(f"\n  -- Occupied bandwidth --")
            print(f"  Analysis band     : {B99_LOW/1e3:.1f} - {B99_HIGH/1e3:.1f} kHz")
            print(f"  B99 edges         : {fl:.0f} - {fh:.0f} Hz")
            print(f"  B99               : {b99:.0f} Hz")
            if kind == "pm":
                carson = 2.0 * (cfg["kp"] + 1.0) * FM
                dev = 100.0 * (b99 - carson) / carson
                flag = "" if abs(dev) <= 20 else "   <-- OUTSIDE 20% TOLERANCE"
                print(f"  Carson estimate   : {carson:.0f} Hz  ({dev:+.1f} %){flag}")
        except ValueError as e:
            print(f"  B99 failed: {e}")

        # ---------- shoulders ----------
        p_up = lobe_power(f, P, FC + SHOULDER_OFFSET)
        p_dn = lobe_power(f, P, FC - SHOULDER_OFFSET)
        print(f"\n  -- Shoulder attenuation (ref: carrier) --")
        print(f"  +2 kHz (7000 Hz)  : {to_db(p_carrier)-to_db(p_up):8.2f} dB")
        print(f"  -2 kHz (3000 Hz)  : {to_db(p_carrier)-to_db(p_dn):8.2f} dB")
        if kind == "pm" and abs(cfg["kp"] - 2.4048) < 0.01:
            print("  NOTE: carrier is nulled here, so shoulder attenuation")
            print("        referenced to the carrier is not meaningful.")

    # ---------- BPSK ----------
    if kind == "bpsk":
        print(f"\n  -- BPSK main lobe --")
        f_nl, p_nl = find_null(f, P, FC - RS)
        f_nu, p_nu = find_null(f, P, FC + RS)
        print(f"  Lower null        : {f_nl:.0f} Hz   (expected {FC-RS:.0f})")
        print(f"  Upper null        : {f_nu:.0f} Hz   (expected {FC+RS:.0f})")
        if np.isfinite(f_nl) and np.isfinite(f_nu):
            print(f"  Null-to-null      : {f_nu-f_nl:.0f} Hz   "
                  f"(expected {2*RS:.0f})")

        f_pk, p_pk = find_peak(f, P, FC - RS, FC + RS)
        f_s1u, p_s1u = find_peak(f, P, FC + RS, FC + 2 * RS)
        f_s1l, p_s1l = find_peak(f, P, FC - 2 * RS, FC - RS)
        print(f"\n  -- Sidelobes (theory: -13.3 dBc for rectangular NRZ) --")
        print(f"  Main lobe peak    : {f_pk:.0f} Hz  {to_db(p_pk):8.2f} dB(V^2)")
        print(f"  Upper 1st sidelobe: {f_s1u:.0f} Hz  "
              f"{to_db(p_s1u)-to_db(p_pk):+8.2f} dBc")
        print(f"  Lower 1st sidelobe: {f_s1l:.0f} Hz  "
              f"{to_db(p_s1l)-to_db(p_pk):+8.2f} dBc")

        print("\n  NOTE: B99% is not reported for BPSK. True B99 for")
        print("        rectangular NRZ is about 20*Rs = 10 kHz, far outside")
        print("        the analysis band. Report null-to-null width and")
        print("        sidelobe level instead.")

    # ---------- floor check ----------
    if floor is not None:
        med_sig = to_db(np.median(P))
        med_flr = to_db(np.median(floor["Pxx"]))
        print(f"\n  -- Noise floor --")
        print(f"  Median PSD, signal: {med_sig:8.2f} dB(V^2/Hz)")
        print(f"  Median PSD, floor : {med_flr:8.2f} dB(V^2/Hz)")
        print(f"  Margin            : {med_sig-med_flr:8.2f} dB")


# ======================================================================
# PLOTS
# ======================================================================

def plot_case(name, res, floor):
    f, P = res["f"], res["Pxx"]
    cfg = res["cfg"]
    fam = cfg["family"]
    Pdb = to_db(P)

    plt.figure(figsize=(9, 5))
    plt.plot(f / 1e3, Pdb, linewidth=0.8, label=f"{name} (Point {cfg['point']})")
    if floor is not None and len(floor["f"]) == len(f):
        plt.plot(f / 1e3, to_db(floor["Pxx"]), linewidth=0.7,
                 alpha=0.6, label="instrument floor")

    if fam == "W":
        plt.axvline(95, linestyle="--", alpha=0.4)
        plt.axvline(105, linestyle="--", alpha=0.4)
        plt.xlim(0, FAMILY_FS["W"] / 2e3)
        title = f"{name} - wide span"
    else:
        plt.axvspan(B99_LOW / 1e3, B99_HIGH / 1e3, alpha=0.12,
                    label="B99 analysis band")
        plt.xlim(0, 25)
        title = f"{name} - signal band"

    plt.xlabel("Frequency (kHz)")
    plt.ylabel("PSD (dB(V$^2$/Hz))")
    plt.title(f"{title}  (average of {res['n_used']} records)")
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=8)
    plt.tight_layout()

    if SAVE_PLOTS:
        plt.savefig(f"{name}_psd.png", dpi=PLOT_DPI)
        plt.savefig(f"{name}_psd.pdf")
    plt.close()


def plot_before_after(res_a, res_b, floor):
    """The before/after LPF overlay - Figure 5 of the thesis."""
    if res_a is None or res_b is None:
        return
    if len(res_a["f"]) != len(res_b["f"]):
        print("\nBefore/after plot skipped: different record lengths.")
        return

    f = res_a["f"]
    plt.figure(figsize=(9, 5))
    plt.plot(f / 1e3, to_db(res_a["Pxx"]), linewidth=0.8,
             label="Before LPF (Point A)")
    plt.plot(f / 1e3, to_db(res_b["Pxx"]), linewidth=0.8,
             label="After LPF (Point B)")
    if floor is not None and len(floor["f"]) == len(f):
        plt.plot(f / 1e3, to_db(floor["Pxx"]), linewidth=0.7,
                 alpha=0.6, linestyle=":", label="Instrument noise floor")
    plt.axvline(95, linestyle="--", alpha=0.4)
    plt.axvline(105, linestyle="--", alpha=0.4)
    plt.xlabel("Frequency (kHz)")
    plt.ylabel("PSD (dB(V$^2$/Hz))")
    plt.title("Spectral comparison before and after Sallen-Key LPF")
    plt.xlim(0, FAMILY_FS["W"] / 2e3)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=8)
    plt.tight_layout()
    if SAVE_PLOTS:
        plt.savefig("FIG_before_after_LPF.png", dpi=PLOT_DPI)
        plt.savefig("FIG_before_after_LPF.pdf")
    plt.close()
    print("\nSaved FIG_before_after_LPF.(png|pdf)")


# ======================================================================
# MAIN
# ======================================================================

def main():
    print("=" * 72)
    print("SPECTRAL ANALYSIS")
    print("=" * 72)

    loaded = {}
    for name, cfg in CASES.items():
        res = collect_case(name, cfg)
        if res is None:
            print(f"\n[{name}] no usable files matching '{cfg['pattern']}'")
            continue
        loaded[name] = res

    if not loaded:
        print("\nNothing loaded. Check DATA_DIR and the CASES patterns.")
        sys.exit(1)

    # ---- file screening summary ----
    print("\n" + "-" * 72)
    print("FILE SCREENING")
    print("-" * 72)
    print(f"{'Case':<7}{'Family':<8}{'Point':<7}{'Found':<8}"
          f"{'Dupes':<8}{'Bad':<7}{'Used':<7}{'Bin (Hz)'}")
    for name, r in loaded.items():
        df = r["f"][1] - r["f"][0]
        print(f"{name:<7}{r['cfg']['family']:<8}{r['cfg']['point']:<7}"
              f"{r['n_files']:<8}{r['n_dup']:<8}{r['n_bad']:<7}"
              f"{r['n_used']:<7}{df:.2f}")

    any_dupes = any(r["n_dup"] for r in loaded.values())
    if any_dupes:
        print("\n  ** Duplicates found. Report the USED count in the thesis,")
        print("     not the number of files on the pendrive. **")

    # ---- per-case analysis ----
    for name, r in loaded.items():
        if r["cfg"]["kind"] == "floor":
            continue

        fam = r["cfg"]["family"]
        floor = loaded.get(FLOOR_FOR_FAMILY.get(fam))

        print("\n" + "=" * 72)
        print(f"{name}   family {fam}   Point {r['cfg']['point']}   "
              f"{r['cfg']['kind']}"
              + (f"   kp={r['cfg']['kp']}" if r["cfg"]["kp"] else ""))
        print("=" * 72)

        if fam == "W":
            analyse_family_w(name, r, floor)
        else:
            analyse_family_n(name, r, floor)

        plot_case(name, r, floor)

    # ---- before/after overlay ----
    plot_before_after(loaded.get("S1"), loaded.get("S2"), loaded.get("S8W"))

    print("\n" + "=" * 72)
    print("REMINDERS FOR THE WRITE-UP")
    print("=" * 72)
    print("- State the B99 analysis band (2-8 kHz) explicitly.")
    print("- Report image suppression as a lower bound where the image")
    print("  sits within 3 dB of the instrument floor.")
    print("- Report the number of records actually used, after dedup.")
    print("- Do not report SINAD or ENOB for BPSK cases.")
    print("- The image suppression figure already appears in the SNR")
    print("  chapter. Report the number once; use this chapter for the")
    print("  PSD figures.")
    print("\nPlots written as <case>_psd.(png|pdf)")


if __name__ == "__main__":
    main()
