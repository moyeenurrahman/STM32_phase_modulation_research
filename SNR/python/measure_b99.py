"""
measure_b99.py  -  measure the 99% power bandwidth (B99%) of the N6 BPSK
capture, directly from the CSVs you already have.

Why this exists: the spectral analysis document did not measure B99% for
BPSK, and the SNR document needs it to set the in-band signal mask for
case N6 (BPSK_HALFBW in snr_method_a.py). Rather than guess a number,
this computes B99% from the same averaged PSD that snr_method_a.py uses,
so both documents end up citing the identical measured value - which is
exactly what the methodology document requires ("do not choose a new
band, or the two documents will disagree").

Definition used (standard, and the one implied by "B99%" in the
literature):

    Take the power spectral density. Centre a band on the CARRIER
    frequency (not on the spectral centroid) and grow it symmetrically
    outward, bin by bin, until the power enclosed reaches 99% of the
    total power in the search window. The full width of that band is
    B99%; half of it is what SIGNAL_MODE = "bpsk" calls BPSK_HALFBW.

    Centring on the carrier, not the centroid, matters here because a
    symmetric digital baseband spectrum should be centred on the carrier
    by construction; centring on the carrier also matches how the SNR
    script places its signal mask (FC +/- BPSK_HALFBW), so the band this
    script measures is the band that will actually be used.

Usage:
    1. Set FOLDER to your N6 folder.
    2. Set SEARCH_HALFBW to a width comfortably larger than you expect
       B99% to be, so the search has room - the default of 20 kHz is
       generous for 500 bps NRZ BPSK.
    3. Run.
    4. Copy the printed BPSK_HALFBW line into snr_method_a.py.
"""

import os
import glob

import numpy as np
import matplotlib.pyplot as plt

# =============================================================
# CONFIGURATION
# =============================================================

FOLDER = "snr_data/N6"
FC     = 5000.0          # carrier frequency, Hz

SEARCH_HALFBW = 40000.0  # how far out to search for the 99% point, Hz.
                          # Must be well inside the capture's Nyquist, and
                          # wide enough that the answer has converged (see
                          # the convergence check the script runs below).
                          # NRZ BPSK has a slowly-decaying sinc^2 spectrum,
                          # so this needs to be much wider than B99% itself.
                          # 40 kHz is confirmed to converge for 500 bps NRZ.

TARGET_FRACTIONS = (0.90, 0.95, 0.99, 0.999)   # printed for context

# Excluded from the power total before computing B99%, same reasoning as
# snr_method_a.py: these are not part of the modulated signal and would
# bias the band outward if left in.
DC_LIMIT = 100.0                    # Hz, DC and very-low-frequency excluded
IMAGES   = (95000.0, 105000.0)      # DAC images, excluded with this half-width
IMG_HALFBW = 1000.0                 # Hz

FS_FALLBACK = 312.5e3     # only used if a CSV has no time column

SAVE_PLOTS = True
PLOT_DPI   = 300


# =============================================================
# LOADER  (same rules as snr_method_a.py, de-duplicated glob)
# =============================================================

def load_csv(path, fs_fallback=None):
    two_col, one_col = [], []
    with open(path, "r", errors="ignore") as fh:
        for line in fh:
            nums, n_fields = [], 0
            for field in line.strip().split(","):
                field = field.strip()
                if not field:
                    continue
                n_fields += 1
                try:
                    nums.append(float(field))
                except ValueError:
                    continue
            if len(nums) >= 2:
                two_col.append((nums[-2], nums[-1]))
            elif len(nums) == 1 and n_fields == 1:
                one_col.append(nums[0])

    if len(two_col) >= 100:
        arr = np.asarray(two_col, dtype=float)
        t, v = arr[:, 0], arr[:, 1]
        if np.median(np.diff(t)) <= 0:
            t, v = v, t
        order = np.argsort(t)
        return t[order], v[order]

    if len(one_col) >= 100:
        if fs_fallback is None:
            raise ValueError(f"{path}: no time column, set FS_FALLBACK")
        v = np.asarray(one_col, dtype=float)
        t = np.arange(len(v)) / float(fs_fallback)
        return t, v

    raise ValueError(f"{path}: could not find data rows")


def _blackmanharris(n):
    a = [0.35875, 0.48829, 0.14128, 0.01168]
    k = np.arange(n)
    return (a[0]
            - a[1] * np.cos(2 * np.pi * k / (n - 1))
            + a[2] * np.cos(4 * np.pi * k / (n - 1))
            - a[3] * np.cos(6 * np.pi * k / (n - 1)))


def psd_one_record(v, fs):
    n = len(v)
    w = _blackmanharris(n)
    vw = (v - np.mean(v)) * w
    spec = np.fft.rfft(vw)
    scale = 2.0 / (fs * np.sum(w ** 2))
    pxx = (np.abs(spec) ** 2) * scale
    pxx[0] /= 2.0
    if n % 2 == 0:
        pxx[-1] /= 2.0
    f = np.fft.rfftfreq(n, 1.0 / fs)
    return f, pxx


def averaged_psd(folder):
    files = glob.glob(os.path.join(folder, "*.CSV")) + \
            glob.glob(os.path.join(folder, "*.csv"))
    seen, uniq = set(), []
    for p in sorted(files):
        key = os.path.normcase(os.path.realpath(p))
        if key not in seen:
            seen.add(key)
            uniq.append(p)
    files = uniq
    if not files:
        raise FileNotFoundError(f"No CSV files found in {folder}")

    acc, f_ref, fs_list = None, None, []
    for path in files:
        t, v = load_csv(path, FS_FALLBACK)
        dt = np.median(np.diff(t))
        fs = 1.0 / dt
        f, pxx = psd_one_record(v, fs)
        if acc is None:
            acc, f_ref = pxx, f
        else:
            if len(pxx) != len(acc):
                raise ValueError(f"{path}: record length differs")
            acc = acc + pxx
        fs_list.append(fs)

    return f_ref, acc / len(files), float(np.mean(fs_list)), len(files)


# =============================================================
# RUN
# =============================================================

print("=" * 62)
print(f"B99% power bandwidth  |  {FOLDER}")
print("=" * 62)

f, pxx, fs, n_rec = averaged_psd(FOLDER)
df = f[1] - f[0]
nyquist = fs / 2.0

print(f"Records averaged : {n_rec}")
print(f"Sample rate      : {fs/1e3:.2f} kSa/s")
print(f"Bin width        : {df:.2f} Hz")
print(f"Nyquist          : {nyquist/1e3:.1f} kHz")

if SEARCH_HALFBW > 0.9 * nyquist:
    raise ValueError("SEARCH_HALFBW is too close to Nyquist - reduce it "
                      "or use a slower timebase.")

# ---------- exclusion masks ----------
dc_mask = f < DC_LIMIT
img_mask = np.zeros_like(f, dtype=bool)
for f_img in IMAGES:
    if f_img < nyquist:
        img_mask |= np.abs(f - f_img) <= IMG_HALFBW

usable = ~dc_mask & ~img_mask

# ---------- power vs half-width, centred on the carrier ----------
search_mask = (np.abs(f - FC) <= SEARCH_HALFBW) & usable
if not np.any(search_mask):
    raise RuntimeError("No usable bins in the search window - check FC "
                       "and SEARCH_HALFBW.")

f_search   = f[search_mask]
pxx_search = pxx[search_mask]
half_widths = np.abs(f_search - FC)

order = np.argsort(half_widths)
half_widths_sorted = half_widths[order]
power_sorted        = pxx_search[order] * df       # power in each bin

cum_power  = np.cumsum(power_sorted)
total_power = cum_power[-1]
cum_frac   = cum_power / total_power

print(f"\nTotal power in search window (+/-{SEARCH_HALFBW/1e3:.0f} kHz "
      f"around carrier, DC and images excluded): "
      f"{total_power:.4e} V^2")

# ---------- find the half-width for each target fraction ----------
#
# IMPORTANT: total_power above is defined as the power WITHIN this search
# window. That means cum_frac always reaches exactly 100% at the window's
# own edge, by construction - so checking "did cum_frac reach 99%?" can
# never fail and proves nothing about whether SEARCH_HALFBW was wide
# enough. The real question is whether the answer is STABLE if the window
# is widened. That is checked below by rerunning with 2x the window.

def b99_halfwidth_for(half_bw):
    """Compute the 99% half-width using a search window of +/-half_bw."""
    m = (np.abs(f - FC) <= half_bw) & usable
    if not np.any(m):
        return None, 0.0
    hw_all = np.abs(f[m] - FC)
    p_all  = pxx[m] * df
    o = np.argsort(hw_all)
    hw_s, p_s = hw_all[o], p_all[o]
    cf = np.cumsum(p_s) / np.sum(p_s)
    idx = np.searchsorted(cf, 0.99)
    idx = min(idx, len(hw_s) - 1)
    return hw_s[idx], hw_s[-1]

print("\nCumulative power vs half-width from the carrier "
      f"(window = +/-{SEARCH_HALFBW/1e3:.1f} kHz):")
for frac in TARGET_FRACTIONS:
    idx = np.searchsorted(cum_frac, frac)
    idx = min(idx, len(half_widths_sorted) - 1)
    hw = half_widths_sorted[idx]
    print(f"  {frac*100:5.1f}% of power within +/-{hw/1e3:6.3f} kHz "
          f"of the carrier  (full bandwidth {2*hw/1e3:6.3f} kHz)")

idx99 = np.searchsorted(cum_frac, 0.99)
idx99 = min(idx99, len(half_widths_sorted) - 1)
b99_halfwidth = half_widths_sorted[idx99]
b99_full = 2 * b99_halfwidth

# ---------- convergence check: does doubling the window change the answer? ----------
check_bw = min(2.0 * SEARCH_HALFBW, 0.95 * nyquist)
if check_bw > SEARCH_HALFBW * 1.05:
    hw_check, edge_reached = b99_halfwidth_for(check_bw)
    if hw_check is not None:
        rel_change = abs(hw_check - b99_halfwidth) / max(b99_halfwidth, 1e-9)
        print(f"\nConvergence check: repeated the search at "
              f"+/-{check_bw/1e3:.1f} kHz (double the window).")
        print(f"  B99% half-width at this window : {hw_check/1e3:.3f} kHz")
        print(f"  Change from the original result : {rel_change*100:.1f}%")
        if rel_change > 0.02:
            print()
            print("  *** NOT CONVERGED ***")
            print(f"  Doubling SEARCH_HALFBW moved the B99% answer by "
                  f"{rel_change*100:.1f}%.")
            print("  This means SEARCH_HALFBW is too narrow to trust the")
            print("  result - there is real signal power beyond the")
            print("  original window that is being missed.")
            print(f"  Increase SEARCH_HALFBW well beyond "
                  f"{check_bw/1e3:.1f} kHz and rerun.")
            print("=" * 62)
            raise SystemExit(1)
        else:
            print("  Converged - widening the window did not move the "
                  "answer. SEARCH_HALFBW is wide enough.")
else:
    print(f"\nNote: SEARCH_HALFBW is already close to Nyquist "
          f"({nyquist/1e3:.0f} kHz), so the convergence check could not "
          "double it. Treat the result with more caution and consider a "
          "slower timebase if possible.")

# ---------- sanity check against the theoretical estimate ----------
# For rectangular NRZ pulses the commonly cited range for 99% power
# bandwidth is roughly 8-20 x the symbol rate, depending on definition.
# This is only a plausibility check, not a substitute for the measurement.
RS_ASSUMED = 500.0
theory_lo, theory_hi = 8 * RS_ASSUMED, 20 * RS_ASSUMED
print(f"\nPlausibility check: rectangular NRZ at {RS_ASSUMED:.0f} bps "
      f"typically gives B99% in the range {theory_lo/1e3:.1f}-"
      f"{theory_hi/1e3:.1f} kHz (rule of thumb, not a strict bound).")
if not (theory_lo <= b99_full <= theory_hi):
    print("  Measured value falls outside that rough range. Not "
          "necessarily wrong - check the plot below - but worth a look.")
else:
    print("  Measured value falls within the expected range.")

# =============================================================
# RESULT
# =============================================================

print("\n" + "-" * 62)
print("RESULT")
print("-" * 62)
print(f"B99% (full bandwidth)     : {b99_full/1e3:.3f} kHz "
      f"({b99_full:.0f} Hz)")
print(f"B99% (half-width, BPSK_HALFBW) : {b99_halfwidth:.0f} Hz")
print()
print("Paste this into snr_method_a.py:")
print()
print(f'    SIGNAL_MODE = "bpsk"')
print(f'    BPSK_HALFBW = {b99_halfwidth:.1f}')
print("=" * 62)

# =============================================================
# PLOTS
# =============================================================

fig, ax = plt.subplots(2, 1, figsize=(10, 8))

pxx_dB = 10 * np.log10(pxx + 1e-30)
ax[0].plot(f / 1e3, pxx_dB, linewidth=0.7)
ax[0].axvspan((FC - b99_halfwidth) / 1e3, (FC + b99_halfwidth) / 1e3,
              alpha=0.25, color="tab:green", label="B99% band")
ax[0].axvline(FC / 1e3, linestyle="--", alpha=0.4, color="k")
ax[0].set_xlim(max(0, (FC - SEARCH_HALFBW) / 1e3),
               (FC + SEARCH_HALFBW) / 1e3)
ax[0].set_xlabel("Frequency (kHz)")
ax[0].set_ylabel("PSD (dB/Hz)")
ax[0].set_title(f"N6 BPSK spectrum - B99% = {b99_full/1e3:.3f} kHz "
                f"(dashed = carrier)")
ax[0].legend(fontsize=9)
ax[0].grid(True, alpha=0.3)

ax[1].plot(half_widths_sorted / 1e3, cum_frac * 100, linewidth=1.2)
ax[1].axhline(99.0, linestyle="--", color="tab:red", alpha=0.6,
              label="99%")
ax[1].axvline(b99_halfwidth / 1e3, linestyle="--", color="tab:green",
              alpha=0.6, label=f"+/-{b99_halfwidth/1e3:.3f} kHz")
ax[1].set_xlabel("Half-width from carrier (kHz)")
ax[1].set_ylabel("Cumulative power (%)")
ax[1].set_title("Cumulative power vs distance from the carrier")
ax[1].set_ylim(0, 101)
ax[1].legend(fontsize=9)
ax[1].grid(True, alpha=0.3)

plt.tight_layout()
if SAVE_PLOTS:
    out = f"b99_{os.path.basename(FOLDER)}.png"
    plt.savefig(out, dpi=PLOT_DPI)
    print(f"Saved {out}")
plt.show()
