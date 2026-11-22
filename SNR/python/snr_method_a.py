"""
Method A - Frequency-domain SNR, SFDR and DAC image levels.

Written for the Tektronix TBS1000C, which gives 2000-point records.

The methodology document assumes 262144-point records and uses Welch's method
with nperseg = 65536 inside a single capture. That is impossible here. Instead,
this script averages the power spectral density across many separate short
records. The statistical effect is the same: averaging M records reduces the
variance of the noise estimate by the same factor as M Welch segments would.

Usage:
    1. Put the CSV files for one condition in one folder.
    2. Set FOLDER and FLOOR_FOLDER below.
    3. Run.
    4. Record the printed numbers in Table 1.

Run once per condition (N2 point A, N3 point B, and so on).
"""

import os
import glob

import numpy as np
import matplotlib.pyplot as plt

# =============================================================
# CONFIGURATION - edit these for each condition
# =============================================================

FOLDER       = "snr_data/N3_pointB_S1"   # condition being analysed
FLOOR_FOLDER = "snr_data/N1_floor_S1"    # instrument noise floor, case N1
LABEL        = "N3 - Point B (LPF output)"

FC = 5000.0              # carrier frequency, Hz

# ---------- which kind of signal is in this folder ----------
#
#   "tone"  - unmodulated carrier.  Cases N1, N2, N3.
#             Signal band is +/-SIG_BINS around the carrier.
#             This is the IEEE 1241 SNR, and it is the headline result.
#
#   "pm"    - analog phase modulation.  Cases N4, N5.
#             Signal band is set by Carson's rule, using KP below.
#
#   "bpsk"  - BPSK.  Case N6.
#             Signal band is +/-BPSK_HALFBW, which should come from the
#             measured B99% in the spectral analysis document.
#
# IMPORTANT: for "pm" and "bpsk" the result is NOT IEEE 1241 SNR. The
# sidebands of a phase-modulated signal are given by Bessel functions and
# extend outwards indefinitely, so any signal/noise boundary is a choice
# rather than a definition. These cases are reported as in-band SNR with
# the band stated explicitly, and are useful only for comparing operating
# conditions against each other.
SIGNAL_MODE = "tone"

# --- for SIGNAL_MODE = "pm" ---
KP = 1.0                 # phase deviation in radians. N4: 1.0, N5: 2.4048
FM = 500.0               # modulating frequency, Hz

# --- for SIGNAL_MODE = "bpsk" ---
BPSK_HALFBW = 1500.0     # half-width of the signal band, Hz.
                         # Replace with the measured B99%/2 from the
                         # spectral analysis document. The default here is
                         # 3 x the 500 bps symbol rate, which is only a
                         # placeholder.

DC_LIMIT   = 100.0       # everything below this is excluded
INBAND_TOP = 50000.0     # in-band noise limit = DAC Nyquist
WIDE_TOP   = 200000.0    # wideband limit; reduced automatically if the
                         # capture's Nyquist is lower

# Band half-widths expressed in BINS, not Hz.
# A 4-term Blackman-Harris window has a main lobe 8 bins wide, so +/-4 bins
# captures the whole tone. 5 is used for margin.
SIG_BINS  = 5            # signal band half-width, "tone" mode only
HARM_BINS = 5            # harmonic exclusion half-width
IMG_BINS  = 6            # image exclusion half-width

NHARM  = 10              # exclude harmonics 2..10 (see note below)
IMAGES = (95000.0, 105000.0)   # zero-order-hold images at fs_dac +/- fc

FLOOR_MARGIN_DB = 3.0    # how close to the floor counts as instrument limited

# Only used if the CSV has no time column (some OpenChoice exports).
# Must match the sample rate shown on the oscilloscope.
FS_FALLBACK = 312.5e3

SAVE_PLOTS = True
PLOT_DPI   = 300

# Note on harmonic exclusion.
# The carrier sits at exactly fs_dac/20, so harmonics above 50 kHz alias back
# onto multiples of 5 kHz rather than landing at new frequencies. Excluding
# only harmonics 2-5 would leave aliased energy counted as noise and would
# depress the reported SNR. The loop must run from 2 to 10.


# =============================================================
# CSV LOADER
# =============================================================

def load_tek_csv(path, fs_fallback=None):
    """
    Read a CSV and return (time, voltage) as float arrays.

    Handles the three layouts you can get from a Tektronix scope:

      1. Classic Tek pendrive format
             Record Length,2.0E3,,-4.0E-03,-1.60E-02,
         (metadata in the first columns, data in the last two)

      2. OpenChoice Desktop two-column export
             TIME,CH1
             -4.000E-03,-1.600E-02

      3. Voltage only, no time column
         Some OpenChoice settings export just the samples. In that case
         the time axis is rebuilt from FS_FALLBACK, so that value must
         match the sample rate shown on the oscilloscope.

    The layout is detected automatically. Metadata lines are ignored
    because they contain fewer than two numbers.
    """
    two_col, one_col = [], []

    with open(path, "r", errors="ignore") as fh:
        for line in fh:
            nums = []
            n_fields = 0
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
                # Accept a bare number on its own line only. A line like
                # "Record Length,2000" also contains exactly one number,
                # but it is metadata, not a sample - the extra text field
                # is what tells the two apart.
                one_col.append(nums[0])

    if len(two_col) >= 100:
        arr = np.asarray(two_col, dtype=float)
        t, v = arr[:, 0], arr[:, 1]

        # If the first column is not increasing it is not a time axis.
        # Some exports put the voltage first.
        if np.median(np.diff(t)) <= 0:
            t, v = v, t
        if np.median(np.diff(t)) <= 0:
            raise ValueError(f"{path}: neither column looks like a time axis")

        order = np.argsort(t)
        return t[order], v[order]

    if len(one_col) >= 100:
        if fs_fallback is None:
            raise ValueError(
                f"{path}: no time column found. Set FS_FALLBACK to the "
                "sample rate shown on the oscilloscope and run again."
            )
        v = np.asarray(one_col, dtype=float)
        t = np.arange(len(v), dtype=float) / float(fs_fallback)
        return t, v

    raise ValueError(
        f"{path}: only {len(two_col)} two-column and {len(one_col)} "
        "single-column data rows found. The export was probably empty."
    )


def psd_one_record(v, fs):
    """
    One-sided power spectral density of a single record, using a
    Blackman-Harris window and density scaling.

    Blackman-Harris is used rather than Hann because the STM32 and the
    oscilloscope run from independent clocks, so the carrier never lands
    exactly on a bin. Hann side lobes near -31 dB would leak carrier energy
    into the noise band and depress the measured SNR badly. Blackman-Harris
    side lobes are near -92 dB.
    """
    n = len(v)
    w = _blackmanharris(n)

    vw = (v - np.mean(v)) * w
    spec = np.fft.rfft(vw)

    # density scaling, one-sided
    scale = 2.0 / (fs * np.sum(w ** 2))
    pxx = (np.abs(spec) ** 2) * scale
    pxx[0] /= 2.0
    if n % 2 == 0:
        pxx[-1] /= 2.0

    f = np.fft.rfftfreq(n, 1.0 / fs)
    return f, pxx


def _blackmanharris(n):
    """4-term Blackman-Harris window."""
    a = [0.35875, 0.48829, 0.14128, 0.01168]
    k = np.arange(n)
    return (a[0]
            - a[1] * np.cos(2 * np.pi * k / (n - 1))
            + a[2] * np.cos(4 * np.pi * k / (n - 1))
            - a[3] * np.cos(6 * np.pi * k / (n - 1)))


def averaged_psd(folder):
    """
    Average the PSD across every CSV in a folder.
    Returns (f, mean_psd, fs, n_records, n_points).
    """
    # Case-insensitive filesystems (Windows) match BOTH patterns for the
    # same file, so de-duplicate by real path or every record is counted
    # twice.
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

    acc = None
    fs_list, n_list = [], []

    for path in files:
        t, v = load_tek_csv(path, FS_FALLBACK)
        dt = np.median(np.diff(t))
        fs = 1.0 / dt
        f, pxx = psd_one_record(v, fs)

        if acc is None:
            acc = pxx
            f_ref = f
        else:
            if len(pxx) != len(acc):
                raise ValueError(
                    f"{path} has a different record length. "
                    "All records in a folder must use identical settings."
                )
            acc = acc + pxx

        fs_list.append(fs)
        n_list.append(len(v))

    mean_psd = acc / len(files)
    fs_mean = float(np.mean(fs_list))

    spread = (max(fs_list) - min(fs_list)) / fs_mean
    if spread > 0.01:
        print("  WARNING: sample rate varies by more than 1% across records.")
        print("           The timebase was probably changed mid-session.")

    return f_ref, mean_psd, fs_mean, len(files), int(np.median(n_list))


# =============================================================
# LOAD
# =============================================================

print("=" * 62)
print(f"SNR - Method A   |   {LABEL}")
print("=" * 62)

f, pxx, fs, n_rec, n_pts = averaged_psd(FOLDER)
df = f[1] - f[0]
nyquist = fs / 2.0

print(f"Records averaged   : {n_rec}")
print(f"Points per record  : {n_pts}")
print(f"Sample rate        : {fs/1e3:.2f} kSa/s")
print(f"Record duration    : {n_pts/fs*1e3:.3f} ms")
print(f"Bin width          : {df:.2f} Hz")
print(f"Nyquist            : {nyquist/1e3:.1f} kHz")

# --- sanity checks -------------------------------------------------
# Work out where each DAC image actually lands after aliasing. On a scope
# with a fixed 2000-point record the sample rate is set by the timebase, and
# some timebases fold an image straight into the measurement band.
#
# On the TBS1000C:  fs = 125 / (timebase in ms) kSa/s
#   400 us/div -> 312.5 kSa/s : images stay at 95 and 105 kHz.  Correct.
#   1 ms/div   -> 125   kSa/s : 95 kHz folds to 30 kHz, inside the noise
#                               band at -25.6 dBc. Unusable for SNR.
for f_img in IMAGES:
    folded = abs(((f_img + fs / 2.0) % fs) - fs / 2.0)
    if abs(folded - f_img) > 1.0:            # it aliased
        print()
        if folded < INBAND_TOP:
            print("  *** STOP ***")
            print(f"  The {f_img/1e3:.0f} kHz DAC image folds to "
                  f"{folded/1e3:.1f} kHz, inside the noise band.")
            print("  The in-band SNR from this capture is not valid.")
            print("  Recapture at 400 us/div (312.5 kSa/s).")
        else:
            print(f"  Note: the {f_img/1e3:.0f} kHz image folds to "
                  f"{folded/1e3:.1f} kHz.")
            print("  Outside the in-band region, but image levels read from")
            print("  this capture are not trustworthy.")
        print()

wide_top = min(WIDE_TOP, 0.95 * nyquist)
if wide_top < WIDE_TOP:
    print(f"  Note: wideband limit reduced to {wide_top/1e3:.0f} kHz "
          f"(Nyquist is {nyquist/1e3:.0f} kHz).")

inband_top = min(INBAND_TOP, 0.95 * nyquist)


# =============================================================
# BAND MASKS
# =============================================================

harm_hw = HARM_BINS * df
img_hw  = IMG_BINS  * df

mode = SIGNAL_MODE.lower()
if mode not in ("tone", "pm", "bpsk"):
    raise ValueError(f"SIGNAL_MODE must be tone, pm or bpsk - got {SIGNAL_MODE}")

# ---------- signal band ----------
if mode == "tone":
    sig_hw = SIG_BINS * df
    band_note = f"+/-{SIG_BINS} bins"
    ieee_valid = True

elif mode == "pm":
    # Carson's rule: B = 2*(kp + 1)*fm, so the half-width is (kp + 1)*fm.
    #
    # The methodology document writes the band as fc +/- (kp + 2)*fm, which
    # is wider than Carson. The wider band is the safer choice here: it is
    # more forgiving of where the sideband energy actually falls, and using
    # a band that is too narrow would push real signal power into the noise
    # term and understate the SNR. The document's form is kept so the two
    # agree, and the strict Carson figure is printed alongside for reference.
    sig_hw = (KP + 2.0) * FM
    carson_hw = (KP + 1.0) * FM
    band_note = f"kp={KP:g}, fm={FM:g} Hz"
    ieee_valid = False

else:  # bpsk
    sig_hw = BPSK_HALFBW
    band_note = "from measured B99%"
    ieee_valid = False

# guard: the band must contain enough bins to mean anything
if sig_hw < 2.0 * df:
    print(f"  WARNING: signal band is only +/-{sig_hw/df:.1f} bins wide.")
    print("           Too narrow to be reliable. Use a slower timebase.")

signal_mask = np.abs(f - FC) <= sig_hw
dc_mask     = f < DC_LIMIT

# ---------- harmonic exclusion ----------
# Only meaningful for an unmodulated tone. For PM and BPSK the energy near
# a "harmonic" is a mixture of real harmonics and modulation products that
# cannot be separated, so nothing is excluded and the harmonics are simply
# counted as noise. Say so in the results rather than hiding it.
harmonic_mask = np.zeros_like(f, dtype=bool)
if mode == "tone":
    for k in range(2, NHARM + 1):
        harmonic_mask |= np.abs(f - k * FC) <= harm_hw

image_mask = np.zeros_like(f, dtype=bool)
for f_img in IMAGES:
    if f_img < nyquist:
        image_mask |= np.abs(f - f_img) <= img_hw

print(f"\nSignal mode        : {mode.upper()}"
      f"{'' if ieee_valid else '   (NOT IEEE 1241 SNR)'}")
print(f"Signal band        : {FC-sig_hw:.0f} - {FC+sig_hw:.0f} Hz "
      f"({band_note})")
if mode == "pm":
    print(f"  strict Carson    : {FC-carson_hw:.0f} - {FC+carson_hw:.0f} Hz "
          f"(B = 2(kp+1)fm = {2*carson_hw:.0f} Hz)")
    print(f"  band used is wider, per the methodology document")
if mode != "tone":
    print("  harmonics NOT excluded - they cannot be separated from")
    print("  modulation products, so they are counted as noise")

if not np.any(signal_mask):
    raise RuntimeError("Signal band contains no bins. Check FC and the data.")

# ---------- locate the peak ----------
i_peak = np.argmax(pxx * signal_mask)
f_peak = f[i_peak]

if mode == "tone":
    print(f"Carrier peak found : {f_peak:.1f} Hz "
          f"(offset {f_peak - FC:+.1f} Hz)")
    if abs(f_peak - FC) > 2 * df:
        print("  WARNING: carrier is not where expected. "
              "Check the firmware build.")
else:
    # For PM the strongest line in the band need not be the carrier. At
    # kp = 2.4048 (case N5) the carrier sits at the first zero of J0, so
    # the peak will normally be a sideband instead. That is the expected
    # result, not an error - report it and move on.
    print(f"Strongest line     : {f_peak:.1f} Hz "
          f"(offset {f_peak - FC:+.1f} Hz from the carrier)")

    i_fc = int(np.argmin(np.abs(f - FC)))
    carrier_rel = 10 * np.log10(pxx[i_fc] / pxx[i_peak] + 1e-30)
    print(f"Carrier bin level  : {carrier_rel:+.1f} dB relative to the "
          "strongest line")
    if carrier_rel < -10.0:
        print("  Carrier is strongly suppressed - consistent with operation")
        print("  at or near a Bessel null (expected for kp = 2.4048).")


# =============================================================
# IN-BAND SNR
# =============================================================
# Many bins are removed by the harmonic and image exclusions, so noise power
# is estimated as (mean noise density) x (full band width) rather than a raw
# sum over the surviving bins. A raw sum would understate the noise in
# proportion to how many bins were excluded.

inband_region = (f >= DC_LIMIT) & (f <= inband_top)
noise_inband  = inband_region & ~signal_mask & ~harmonic_mask & ~image_mask

P_signal = np.sum(pxx[signal_mask]) * df

band_width  = inband_top - DC_LIMIT
noise_dens  = np.mean(pxx[noise_inband])
P_noise_in  = noise_dens * band_width

n_kept  = int(np.sum(noise_inband))
n_total = int(np.sum(inband_region))

SNR_inband_dB = 10 * np.log10(P_signal / P_noise_in)

# robust variant: drop bins more than 10x the median, in case a spur was missed
med = np.median(pxx[noise_inband])
clean = noise_inband & (pxx <= 10 * med)
P_noise_in_robust = np.mean(pxx[clean]) * band_width
SNR_inband_robust = 10 * np.log10(P_signal / P_noise_in_robust)
n_spurs = n_kept - int(np.sum(clean))


# =============================================================
# WIDEBAND SNR
# =============================================================
# Here the images ARE counted as noise - that is the point of this figure.
# Direct summation is used because few bins are excluded.

wide_region = (f >= DC_LIMIT) & (f <= wide_top)
noise_wide  = wide_region & ~signal_mask
P_noise_wide = np.sum(pxx[noise_wide]) * df
SNR_wide_dB  = 10 * np.log10(P_signal / P_noise_wide)


# =============================================================
# SFDR AND IMAGE LEVELS
# =============================================================

pxx_dB = 10 * np.log10(pxx + 1e-30)
carrier_dB = np.max(pxx_dB[signal_mask])

spur_region = wide_region & ~signal_mask & ~dc_mask
sfdr_dB = carrier_dB - np.max(pxx_dB[spur_region])
f_worst = f[spur_region][np.argmax(pxx_dB[spur_region])]


def peak_near(f_target, half_width):
    m = np.abs(f - f_target) <= half_width
    if not np.any(m):
        return None
    return np.max(pxx_dB[m])


image_levels = {}
for f_img in IMAGES:
    if f_img < nyquist:
        pk = peak_near(f_img, img_hw)
        image_levels[f_img] = carrier_dB - pk
    else:
        image_levels[f_img] = None


# =============================================================
# INSTRUMENT FLOOR COMPARISON
# =============================================================

floor_available = os.path.isdir(FLOOR_FOLDER)
verdict = "not checked"
P_noise_floor = None

if floor_available:
    f_fl, pxx_fl, fs_fl, n_fl, _ = averaged_psd(FLOOR_FOLDER)

    if len(f_fl) != len(f) or abs(fs_fl - fs) / fs > 0.01:
        print("\n  WARNING: floor capture used different settings.")
        print("           The comparison is not valid. Recapture the floor.")
        floor_available = False
    else:
        fl_region = (f_fl >= DC_LIMIT) & (f_fl <= inband_top)
        P_noise_floor = np.mean(pxx_fl[fl_region]) * band_width

        excess_dB = 10 * np.log10(P_noise_in / P_noise_floor)

        if excess_dB > FLOOR_MARGIN_DB:
            verdict = "MEASURED VALUE"
        else:
            verdict = "LOWER BOUND - instrument limited"


# =============================================================
# RESULTS
# =============================================================

print()
print("-" * 62)
if ieee_valid:
    print("RESULTS   (IEEE Std 1241 SNR - single-tone test signal)")
else:
    print("RESULTS   *** IN-BAND SNR ONLY - NOT IEEE 1241 SNR ***")
print("-" * 62)
if not ieee_valid:
    print(f"Band used: {FC-sig_hw:.0f} - {FC+sig_hw:.0f} Hz. This value is")
    print("meaningless without that band, so quote them together in the")
    print("thesis table. Use it to compare operating conditions against")
    print("each other, not as an absolute figure of merit.")
    print()
print(f"Signal power           : {P_signal:.4e} V^2  "
      f"({np.sqrt(P_signal)*1e3:.1f} mV rms)")
print(f"In-band noise power    : {P_noise_in:.4e} V^2")
print(f"Noise bins used        : {n_kept} of {n_total} "
      f"({100*n_kept/n_total:.0f}% of the band)")
if n_spurs:
    print(f"Outlier bins dropped   : {n_spurs} (possible missed spurs)")
print()
print(f"In-band SNR  (0.1-{inband_top/1e3:.0f} kHz) : {SNR_inband_dB:7.2f} dB")
print(f"  robust variant                : {SNR_inband_robust:7.2f} dB")
print(f"Wideband SNR (0.1-{wide_top/1e3:.0f} kHz): {SNR_wide_dB:7.2f} dB")
print(f"SFDR                            : {sfdr_dB:7.2f} dB "
      f"(worst spur at {f_worst/1e3:.2f} kHz)")
for f_img, lvl in image_levels.items():
    if lvl is None:
        print(f"Image at {f_img/1e3:.0f} kHz            : above Nyquist, "
              f"not measurable at this timebase")
    else:
        print(f"Image at {f_img/1e3:.0f} kHz            : {-lvl:7.2f} dBc")

print()
print("-" * 62)
print("INSTRUMENT LIMIT CHECK")
print("-" * 62)

if floor_available:
    SNR_floor_limit = 10 * np.log10(P_signal / P_noise_floor)
    print(f"Floor records averaged : {n_fl}")
    print(f"Floor noise power      : {P_noise_floor:.4e} V^2")
    print(f"Measured noise is      : {excess_dB:+.2f} dB above the floor")
    print(f"Scope's own SNR limit  : {SNR_floor_limit:.2f} dB "
          f"(the best this setup could ever show)")
    print()
    print(f"  VERDICT: {verdict}")

    if excess_dB > FLOOR_MARGIN_DB:
        P_corr = P_noise_in - P_noise_floor
        SNR_corr = 10 * np.log10(P_signal / P_corr)
        print(f"  Floor-subtracted estimate: {SNR_corr:.2f} dB")
        print("  (supporting remark only - state the subtraction in the text)")
        print()
        print(f"  Report as: SNR = {SNR_inband_dB:.1f} dB")
    else:
        print("  The oscilloscope noise dominates. The DAC noise is buried.")
        print("  Do NOT report this as a measured value.")
        print()
        print(f"  Report as: SNR >= {SNR_inband_dB:.1f} dB, "
              f"measurement limited by oscilloscope resolution")
else:
    print("  No valid floor capture found.")
    print("  Every SNR value from this run is unverified. Capture case N1.")

print("=" * 62)


# =============================================================
# PLOTS
# =============================================================

fig, ax = plt.subplots(figsize=(11, 6))
ax.plot(f / 1e3, pxx_dB, linewidth=0.7, label=LABEL)

if floor_available:
    ax.plot(f_fl / 1e3, 10 * np.log10(pxx_fl + 1e-30),
            linewidth=0.7, alpha=0.75, label="N1 instrument floor")

ax.axvspan((FC - sig_hw) / 1e3, (FC + sig_hw) / 1e3,
           alpha=0.30, color="tab:green", label="Signal band")
ax.axvspan(DC_LIMIT / 1e3, inband_top / 1e3,
           alpha=0.07, color="tab:blue", label="In-band noise region")

for f_img in IMAGES:
    if f_img < nyquist:
        ax.axvline(f_img / 1e3, linestyle="--", linewidth=0.8,
                   alpha=0.5, color="tab:red")

ax.set_xlabel("Frequency (kHz)")
ax.set_ylabel("PSD (dB/Hz)")
ax.set_title(f"SNR analysis - {LABEL}\n"
             f"{n_rec} records averaged, {df:.0f} Hz bins, "
             f"Blackman-Harris window")
ax.set_xlim(0, wide_top / 1e3)
ax.grid(True, alpha=0.3)
ax.legend(loc="upper right", fontsize=9)
plt.tight_layout()

if SAVE_PLOTS:
    out = f"psd_{os.path.basename(FOLDER)}.png"
    plt.savefig(out, dpi=PLOT_DPI)
    print(f"Saved {out}")

# zoom on the image region - this is the plot that shows the LPF working
if IMAGES[0] < nyquist:
    fig2, ax2 = plt.subplots(figsize=(9, 5))
    ax2.plot(f / 1e3, pxx_dB, linewidth=0.9, label=LABEL)
    if floor_available:
        ax2.plot(f_fl / 1e3, 10 * np.log10(pxx_fl + 1e-30),
                 linewidth=0.9, alpha=0.75, label="N1 instrument floor")
    ax2.set_xlim(80, min(120, nyquist / 1e3))
    ax2.set_xlabel("Frequency (kHz)")
    ax2.set_ylabel("PSD (dB/Hz)")
    ax2.set_title("Zero-order-hold images at 95 kHz and 105 kHz")
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=9)
    plt.tight_layout()
    if SAVE_PLOTS:
        out2 = f"images_{os.path.basename(FOLDER)}.png"
        plt.savefig(out2, dpi=PLOT_DPI)
        print(f"Saved {out2}")

plt.show()
