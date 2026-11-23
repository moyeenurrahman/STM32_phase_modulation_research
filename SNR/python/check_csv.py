"""
check_csv.py  -  inspect ONE oscilloscope CSV before committing to a session.

Answers five questions:
  1. Can the file be parsed at all?
  2. What sample rate did the scope really use?
  3. Is the carrier where it should be?
  4. Is the signal clipping, or too small on the screen?
  5. What does the spectrum look like?

Run this on your first capture. If it passes, the whole pipeline works and
you can go ahead with the full 32-record session.

Put this file in the SNRpython folder, next to the CSV, and run it.
"""

import os
import glob

import numpy as np
import matplotlib.pyplot as plt

# =============================================================
# CONFIGURATION
# =============================================================

FILENAME = "1"          # your file. "1", "1.csv", or a full path all work.
                        # Set to "" to just use the first CSV-ish file found.

FC_EXPECTED = 5000.0    # what the carrier should be, Hz
FS_FALLBACK = 312500.0  # only used if the file has no time column


# =============================================================
# FIND THE FILE
# =============================================================

def find_file(name):
    """Cope with '1', '1.csv', '1.CSV', or no name at all."""
    here = os.path.dirname(os.path.abspath(__file__))

    if name:
        for cand in (name, name + ".csv", name + ".CSV"):
            p = cand if os.path.isabs(cand) else os.path.join(here, cand)
            if os.path.isfile(p):
                return p

    # nothing matched - look for anything that could be a capture
    found = []
    for pat in ("*.csv", "*.CSV", "*"):
        for p in sorted(glob.glob(os.path.join(here, pat))):
            if os.path.isfile(p) and not p.endswith((".py", ".png")):
                found.append(p)
        if found:
            break

    if not found:
        raise FileNotFoundError(
            f"No data file found in {here}\n"
            f"Files present: {os.listdir(here)}"
        )

    print(f"'{name}' not found - using {os.path.basename(found[0])} instead\n")
    return found[0]


# =============================================================
# LOADER  (same rules as snr_method_a.py)
# =============================================================

def load_csv(path, fs_fallback=None):
    """
    Handles the three layouts a Tektronix scope can produce:
      1. Classic pendrive format - metadata columns then time,volts
      2. OpenChoice two-column export
      3. Voltage only, no time column
    """
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
        return t[order], v[order], "time + voltage"

    if len(one_col) >= 100:
        if fs_fallback is None:
            raise ValueError("No time column, and FS_FALLBACK not set.")
        v = np.asarray(one_col, dtype=float)
        t = np.arange(len(v)) / float(fs_fallback)
        return t, v, "voltage only (time rebuilt)"

    raise ValueError(
        f"Could not find data rows. Found {len(two_col)} two-column and "
        f"{len(one_col)} single-column rows."
    )


def blackmanharris(n):
    a = [0.35875, 0.48829, 0.14128, 0.01168]
    k = np.arange(n)
    return (a[0]
            - a[1] * np.cos(2 * np.pi * k / (n - 1))
            + a[2] * np.cos(4 * np.pi * k / (n - 1))
            - a[3] * np.cos(6 * np.pi * k / (n - 1)))


# =============================================================
# RUN
# =============================================================

path = find_file(FILENAME)

print("=" * 60)
print(f"Checking: {os.path.basename(path)}")
print("=" * 60)

# --- show the top of the file, so the format is visible ---
print("\nFirst 6 lines of the file:")
with open(path, "r", errors="ignore") as fh:
    for i, line in enumerate(fh):
        if i >= 6:
            break
        print(f"  {line.rstrip()[:78]}")

t, v, layout = load_csv(path, FS_FALLBACK)

n = len(v)
dt = np.median(np.diff(t))
fs = 1.0 / dt
duration = n * dt
df = 1.0 / duration

print(f"\nLayout detected  : {layout}")
print(f"Points           : {n}")
print(f"Sample rate      : {fs/1e3:.2f} kSa/s")
print(f"Duration         : {duration*1e3:.3f} ms")
print(f"Bin width        : {df:.2f} Hz")
print(f"Nyquist          : {fs/2e3:.1f} kHz")

# --- amplitude ---
vmin, vmax = float(np.min(v)), float(np.max(v))
vpp = vmax - vmin
print(f"\nMinimum          : {vmin:+.4f} V")
print(f"Maximum          : {vmax:+.4f} V")
print(f"Peak to peak     : {vpp:.4f} V")
print(f"Mean (DC)        : {np.mean(v):+.4f} V")

# --- spectrum ---
x = v - np.mean(v)
w = blackmanharris(n)
spec = np.abs(np.fft.rfft(x * w))
freqs = np.fft.rfftfreq(n, dt)

spec_s = spec.copy()
spec_s[freqs < 100.0] = 0.0
f_peak = freqs[np.argmax(spec_s)]

spec_dB = 20 * np.log10(spec + 1e-30)
spec_dBc = spec_dB - np.max(spec_dB)

print(f"\nStrongest peak   : {f_peak:.1f} Hz")
print(f"Expected         : {FC_EXPECTED:.1f} Hz")
print(f"Difference       : {f_peak - FC_EXPECTED:+.1f} Hz "
      f"({(f_peak/FC_EXPECTED - 1)*1e6:+.0f} ppm)")

# --- checks ---
print("\n" + "-" * 60)
print("CHECKS")
print("-" * 60)
ok = True

if n < 500:
    print(f"  FAIL  only {n} points. Export was probably screen data only.")
    ok = False
else:
    print(f"  ok    record length ({n} points)")

# Coherence: for one clean sine captured in a single acquisition, nearly all
# the power sits in the main lobe around the peak. If the scope was still
# running during Get Data, the readout stitches together fragments of many
# acquisitions at random phase, and the power smears across the spectrum.
lobe = np.abs(freqs - f_peak) <= 5 * df
coherence = float(np.sum(spec[lobe] ** 2) / np.sum(spec ** 2))

if coherence < 0.5:
    print(f"  FAIL  only {coherence*100:.0f}% of the power is at one frequency.")
    print("        The record is not a single acquisition. Press Single on the")
    print("        scope and wait for it to stop, THEN press Get Data.")
    ok = False
else:
    print(f"  ok    single coherent acquisition ({coherence*100:.0f}% of power "
          "in the carrier)")

    if abs(f_peak - FC_EXPECTED) > 0.05 * FC_EXPECTED:
        print(f"  FAIL  carrier is {f_peak:.0f} Hz, expected "
              f"~{FC_EXPECTED:.0f} Hz.")
        if abs(f_peak - 2 * FC_EXPECTED) < 0.05 * FC_EXPECTED:
            print("        Exactly double - the DMA data width is set to Word "
                  "instead of Half Word.")
        elif abs(f_peak - 0.5 * FC_EXPECTED) < 0.05 * FC_EXPECTED:
            print("        Exactly half - check the TIM6 clock and ARR value.")
        ok = False
    else:
        print(f"  ok    carrier near {FC_EXPECTED/1e3:.0f} kHz")

# clipping: are many samples sitting at the extremes?
flat_hi = np.sum(v > vmax - 0.002 * vpp)
flat_lo = np.sum(v < vmin + 0.002 * vpp)
if flat_hi > n * 0.005 or flat_lo > n * 0.005:
    print(f"  FAIL  possible clipping ({flat_hi} samples at the top, "
          f"{flat_lo} at the bottom). Do not capture until fixed.")
    ok = False
else:
    print("  ok    peaks look rounded, no clipping")

# aliasing danger
for f_img in (95000.0, 105000.0):
    folded = abs(((f_img + fs / 2) % fs) - fs / 2)
    if folded < 50000.0 and abs(folded - f_img) > 1.0:
        print(f"  WARN  the {f_img/1e3:.0f} kHz DAC image folds to "
              f"{folded/1e3:.1f} kHz, inside the noise band.")
        print(f"        Change the timebase - this corrupts the in-band SNR.")
        ok = False

if ok:
    print("\n  PASS - this file is usable. Go ahead with the full session.")
else:
    print("\n  Fix the problems above before capturing 32 records.")
print("=" * 60)

# =============================================================
# PLOTS
# =============================================================

fig, ax = plt.subplots(2, 1, figsize=(11, 8))

n_show = min(n, int(0.002 / dt))          # first 2 ms
ax[0].plot(t[:n_show] * 1e3, v[:n_show], linewidth=0.8)
ax[0].set_xlabel("Time (ms)")
ax[0].set_ylabel("Voltage (V)")
ax[0].set_title(f"{os.path.basename(path)}  -  first 2 ms  "
                f"({vpp:.3f} Vpp, {f_peak:.0f} Hz)")
ax[0].grid(True, alpha=0.3)

ax[1].plot(freqs / 1e3, spec_dBc, linewidth=0.7)
ax[1].axvline(FC_EXPECTED / 1e3, linestyle="--", alpha=0.5, color="tab:green")
for f_img in (95.0, 105.0):
    if f_img * 1e3 < fs / 2:
        ax[1].axvline(f_img, linestyle="--", alpha=0.4, color="tab:red")
ax[1].set_xlabel("Frequency (kHz)")
ax[1].set_ylabel("dB relative to carrier")
ax[1].set_title(f"Spectrum  -  {df:.0f} Hz bins, Nyquist {fs/2e3:.0f} kHz "
                "(red dashes = DAC images)")
ax[1].set_xlim(0, fs / 2e3)
ax[1].set_ylim(-100, 5)
ax[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "check_csv.png"), dpi=150)
print("Saved check_csv.png")
plt.show()
