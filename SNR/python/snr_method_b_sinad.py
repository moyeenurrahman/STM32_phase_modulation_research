"""
Method B - IEEE Std 1241 four-parameter sine fit: SINAD and ENOB.

Use on the unmodulated carrier only (cases N2 and N3). The sine fit assumes a
single tone. Do not run it on PM or BPSK captures.

Why a four-parameter fit and not a subtraction of an ideal waveform:

  1. The measured waveform has an arbitrary starting phase. Without fitting it,
     the residual can be larger than the signal itself.
  2. The STM32 clock and the oscilloscope clock are independent, so the true
     carrier frequency is never exactly 5000.000 Hz.
  3. The LPF adds -20.7 deg of phase shift at 5 kHz.

Fitting amplitude, frequency, phase and offset absorbs all three problems.

What SINAD measures: the residual after subtracting the fitted sine contains
noise AND harmonic distortion. That quantity is SINAD (equivalently THD+N),
not SNR. Reporting it as SNR would be incorrect.

Important with an 8-bit oscilloscope: the residual will be dominated by the
scope's own quantisation noise, so the ENOB you compute is close to the
OSCILLOSCOPE's effective bits (around 7), not the STM32 DAC's (around 11).
State this in the thesis. It is a limitation of the instrument, not a fault
in the DAC.
"""

import os
import glob

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

# =============================================================
# CONFIGURATION
# =============================================================

FOLDER = "snr_data/N3_pointB_S1"
LABEL  = "N3 - Point B (LPF output)"

FS_DAC    = 100000.0     # DAC sample rate, Hz
DAC_CODES = 2400.0       # peak-to-peak codes used (2 x 1200)
DAC_FULL  = 4096.0       # full scale of a 12-bit DAC

FC_NOMINAL = 5000.0
SEARCH_HZ  = 50.0        # frequency search range around the FFT peak

# Only used if the CSV has no time column (some OpenChoice exports).
# Must match the sample rate shown on the oscilloscope.
FS_FALLBACK = 400e3

SAVE_PLOTS = True
PLOT_DPI   = 300


# =============================================================
# CSV LOADER  (same rule as snr_method_a.py)
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


# =============================================================
# FOUR-PARAMETER FIT
# =============================================================

def fit_record(t, x):
    """
    Fit A*cos(2*pi*f*t) + B*sin(2*pi*f*t) + C to the data.

    For any trial frequency, A, B and C follow from linear least squares.
    So only the frequency needs a search - this is the standard and stable
    way to do the four-parameter fit.
    """
    x = x - np.mean(x)
    dt = np.median(np.diff(t))

    # coarse estimate from the FFT peak
    w = np.blackman(len(x))
    spec = np.abs(np.fft.rfft(x * w))
    freqs = np.fft.rfftfreq(len(x), dt)
    spec[freqs < 100.0] = 0.0
    f0_est = freqs[np.argmax(spec)]

    def design(f0):
        return np.column_stack([
            np.cos(2 * np.pi * f0 * t),
            np.sin(2 * np.pi * f0 * t),
            np.ones_like(t),
        ])

    def fit_at(f0):
        M = design(f0)
        coef, *_ = np.linalg.lstsq(M, x, rcond=None)
        return M @ coef, coef

    def cost(f0):
        fit, _ = fit_at(f0)
        return float(np.sum((x - fit) ** 2))

    res = minimize_scalar(
        cost,
        bounds=(f0_est - SEARCH_HZ, f0_est + SEARCH_HZ),
        method="bounded",
        options={"xatol": 1e-7},
    )

    f0 = float(res.x)
    fit, coef = fit_at(f0)
    A, B, C = coef

    return {
        "f0_coarse": f0_est,
        "f0": f0,
        "amplitude": float(np.hypot(A, B)),
        "phase_deg": float(np.degrees(np.arctan2(-B, A))),
        "offset": float(C),
        "fit": fit,
        "x": x,
        "t": t,
        "residual": x - fit,
    }


# =============================================================
# RUN OVER EVERY RECORD
# =============================================================

# Case-insensitive filesystems (Windows) match BOTH patterns for the same
# file, so de-duplicate by real path or every record is counted twice.
_files = glob.glob(os.path.join(FOLDER, "*.CSV")) + \
         glob.glob(os.path.join(FOLDER, "*.csv"))
_seen, files = set(), []
for _p in sorted(_files):
    _key = os.path.normcase(os.path.realpath(_p))
    if _key not in _seen:
        _seen.add(_key)
        files.append(_p)
if not files:
    raise FileNotFoundError(f"No CSV files found in {FOLDER}")

print("=" * 62)
print(f"SINAD / ENOB - Method B   |   {LABEL}")
print("=" * 62)

correction = 20 * np.log10(DAC_FULL / DAC_CODES)

sinad_list, enob_list, f0_list, amp_list = [], [], [], []
first = None

for path in files:
    t, v = load_tek_csv(path, FS_FALLBACK)
    r = fit_record(t, v)

    P_fit = float(np.mean(r["fit"] ** 2))
    P_res = float(np.mean(r["residual"] ** 2))
    sinad = 10 * np.log10(P_fit / P_res)
    enob = (sinad - 1.76 + correction) / 6.02

    sinad_list.append(sinad)
    enob_list.append(enob)
    f0_list.append(r["f0"])
    amp_list.append(r["amplitude"])

    if first is None:
        first = r

    print(f"  {os.path.basename(path):16s}  "
          f"f0 = {r['f0']:9.3f} Hz   "
          f"A = {r['amplitude']:.4f} V   "
          f"SINAD = {sinad:6.2f} dB")

sinad_list = np.array(sinad_list)
enob_list = np.array(enob_list)
f0_list = np.array(f0_list)
amp_list = np.array(amp_list)

print()
print("-" * 62)
print("RESULTS")
print("-" * 62)
print(f"Records analysed   : {len(files)}")
print(f"Fitted frequency   : {f0_list.mean():.3f} Hz  "
      f"(std {f0_list.std():.3f} Hz)")
print(f"  offset from 5 kHz: {f0_list.mean() - FC_NOMINAL:+.3f} Hz  "
      f"({(f0_list.mean()/FC_NOMINAL - 1)*1e6:+.0f} ppm)")
print(f"Fitted amplitude   : {amp_list.mean():.4f} V  "
      f"(std {amp_list.std()*1e3:.2f} mV)")
print(f"Amplitude correction: {correction:.2f} dB "
      f"({DAC_CODES:.0f} of {DAC_FULL:.0f} codes used)")
print()
print(f"SINAD              : {sinad_list.mean():.2f} dB  "
      f"(std {sinad_list.std():.2f} dB)")
print(f"ENOB               : {enob_list.mean():.2f} bits  "
      f"(std {enob_list.std():.2f})")

print()
print("-" * 62)
print("CHECKS")
print("-" * 62)

# residual structure check
res_first = first["residual"]
t_first = first["t"]
dt = np.median(np.diff(t_first))
w = np.blackman(len(res_first))
spec = np.abs(np.fft.rfft(res_first * w))
freqs = np.fft.rfftfreq(len(res_first), dt)

band = np.abs(freqs - first["f0"]) <= 3 * (freqs[1] - freqs[0])
leak_dB = 20 * np.log10(np.max(spec[band]) / np.median(spec[freqs > 100.0]))

print(f"Residual at the carrier is {leak_dB:.1f} dB above the residual median.")
if leak_dB > 20:
    print("  WARNING: strong structure left at the carrier frequency.")
    print("  The fit did not converge. Discard this result and check the data")
    print("  for clipping or a changed timebase.")
else:
    print("  Residual looks like noise. Fit accepted.")

print()
print("Cross-checks to do by hand:")
print(f"  1. Compare the fitted frequency ({f0_list.mean():.3f} Hz) with your")
print("     output frequency accuracy result. These are independent")
print("     measurements of the same quantity and should agree.")
print("  2. SINAD will come out LOWER than the Method A in-band SNR, for")
print("     two reasons: it includes harmonic distortion, and it is measured")
print("     over the whole captured bandwidth rather than only 0.1-50 kHz.")
print("     At 400 kSa/s that bandwidth difference alone is about 6 dB.")
print("     If SINAD comes out HIGHER, something is wrong in one analysis.")
print("  3. With an 8-bit scope, expect ENOB near 7, not 11. You are")
print("     measuring the instrument, not the DAC.")
print("=" * 62)


# =============================================================
# PLOTS
# =============================================================

n_show = min(len(first["t"]), int(0.002 / dt))   # first 2 ms

fig, ax = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

ax[0].plot(first["t"][:n_show] * 1e3, first["x"][:n_show],
           linewidth=0.9, label="Measured")
ax[0].plot(first["t"][:n_show] * 1e3, first["fit"][:n_show],
           linewidth=0.9, alpha=0.8, label="Fitted sine")
ax[0].set_ylabel("Voltage (V)")
ax[0].set_title(f"Four-parameter sine fit - {LABEL}")
ax[0].legend(fontsize=9)
ax[0].grid(True, alpha=0.3)

ax[1].plot(first["t"][:n_show] * 1e3, first["residual"][:n_show],
           linewidth=0.7, color="tab:red")
ax[1].set_xlabel("Time (ms)")
ax[1].set_ylabel("Residual (V)")
ax[1].set_title("Residual - should look like random noise")
ax[1].grid(True, alpha=0.3)

plt.tight_layout()
if SAVE_PLOTS:
    out = f"sinefit_{os.path.basename(FOLDER)}.png"
    plt.savefig(out, dpi=PLOT_DPI)
    print(f"Saved {out}")

# spectrum of the residual - shows which harmonics remain
fig2, ax2 = plt.subplots(figsize=(10, 5))
ax2.plot(freqs / 1e3, 20 * np.log10(spec + 1e-20), linewidth=0.7)
ax2.set_xlabel("Frequency (kHz)")
ax2.set_ylabel("Residual magnitude (dB)")
ax2.set_title("Spectrum of the residual")
ax2.grid(True, alpha=0.3)
plt.tight_layout()
if SAVE_PLOTS:
    out2 = f"residual_spectrum_{os.path.basename(FOLDER)}.png"
    plt.savefig(out2, dpi=PLOT_DPI)
    print(f"Saved {out2}")

plt.show()
