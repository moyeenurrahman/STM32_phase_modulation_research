"""
make_figures.py   -   thesis figures for the frequency accuracy metric

Generates:
    Fig A   frequency error against test frequency, HSI and HSE
    Fig C   spectrum, unmodulated against beta = 1
    Fig D   error by measurement point (MCO, raw DAC, filter output)
    Fig E   cumulative timing error
    Fig F   uncertainty budget

WHAT IT NEEDS
-------------
Put this file, analyse_freq.py, and a folder of captures together, then
fill in the settings block below. Every figure is optional: if the files
for one are missing, that figure is skipped and the rest still run.

Outputs go to a "figures" sub-folder as PDF (vector, for the thesis) and
PNG (for checking).
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from analyse_freq import load_tek_csv, freq_zero_crossing

# ======================================================================
# SETTINGS
# ======================================================================

OUT = "figures"

# ---- Fig A : from your batch_freq.py grouped output ------------------
# mean and 1 sigma in ppm, in the order 1, 2, 5, 10 kHz
FREQS_KHZ = [1, 2, 5, 10]
HSI_MEAN = [3190.50, 3094.64, 3086.92, 3055.30]
HSI_SD = [114.89, 55.56, 30.56, 96.20]
HSE_MEAN = [65.13, 71.24, 89.19, 105.50]
HSE_SD = [2.42, 7.67, 6.23, 53.03]

# ---- Fig C : one unmodulated and one modulated capture ---------------
FILE_UNMOD = "T2_HSE_5000Hz_run1.csv"
FILE_PM = "T6_HSE_5000Hz_run1.csv"
F_CARRIER = 5000.0
F_MSG = 500.0
BETA = 1.0

# ---- Fig D : HSI results by measurement point ------------------------
# (label, mean ppm, 1 sigma, number of runs)
POINTS = [
    ("MCO pin\n(PA8)",        3047.24, 105.04, 3),
    ("Raw DAC\n(PA4)",        3106.84,  87.15, 12),   # T1 stands in if no
    ("Filter output\n(5 kHz)", 3086.92,  30.56, 3),   # separate raw-DAC set
    ("All frequencies\n(T1)", 3106.84,  87.15, 12),
]
HSI_INSTABILITY = 84.2        # ppm, from analyse_stability.py

# ---- Fig E : one HSI and one HSE capture -----------------------------
FILE_E_HSI = "T1_HSI_5000Hz_run1.csv"
FILE_E_HSE = "T2_HSE_5000Hz_run1.csv"

# ---- Fig F : uncertainty budget --------------------------------------
BUDGET = {
    "HSI": [("Oscillator repeatability", 81.4),
            ("Oscilloscope timebase", 20.0),
            ("Estimator precision", 1.1)],
    "HSE": [("Oscillator repeatability", 17.9),
            ("Oscilloscope timebase", 20.0),
            ("Estimator precision", 1.1)],
}

# ======================================================================

os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({
    "font.size": 9, "axes.grid": True, "grid.alpha": 0.3,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 110,
})
C_HSI, C_HSE, C_REF = "#c0392b", "#2874a6", "#7f8c8d"


def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, name + ".pdf"))
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=200)
    plt.close(fig)
    print(f"  wrote {name}.pdf and .png")


def load(path):
    """Load a capture and return time, volts and the sample interval."""
    d = load_tek_csv(path)
    dt = d["dt"] if d["dt"] else (d["t"][-1] - d["t"][0]) / (len(d["t"]) - 1)
    return d["t"], d["channels"][0], dt


def spectrum(v, dt, ref_hz=None):
    """Blackman-Harris power spectrum in dB relative to the largest line.

    Blackman-Harris is used throughout because the STM32 and the
    oscilloscope run from independent clocks, so the record almost never
    contains a whole number of cycles. Its -92 dB side lobes keep that
    leakage well below the components being measured.
    """
    v = v - np.mean(v)
    N = len(v)
    w = np.blackman(N)                       # numpy's Blackman
    V = np.abs(np.fft.rfft(v * w))
    V[0] = 0.0
    f = np.fft.rfftfreq(N, dt)
    ref = V.max() if ref_hz is None else V[np.argmin(np.abs(f - ref_hz))]
    return f, 20 * np.log10(np.maximum(V, 1e-30) / ref)


# ----------------------------------------------------------------------
# FIGURE A - frequency error against test frequency
# ----------------------------------------------------------------------
print("Fig A")
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0))
x = np.arange(len(FREQS_KHZ))

for a, mean, sd, col, title in (
        (ax[0], HSI_MEAN, HSI_SD, C_HSI, "HSI (internal RC)"),
        (ax[1], HSE_MEAN, HSE_SD, C_HSE, "HSE (crystal reference)")):
    a.errorbar(x, mean, yerr=sd, fmt="o", color=col, capsize=4,
               ms=5, lw=1.2)
    m = np.mean(mean)
    a.axhline(m, color=C_REF, ls="--", lw=1,
              label=f"mean {m:+.0f} ppm")
    a.set_xticks(x)
    a.set_xticklabels([f"{f}" for f in FREQS_KHZ])
    a.set_xlabel("test frequency (kHz)")
    a.set_ylabel("frequency error (ppm)")
    a.set_title(title)
    a.legend(frameon=False, loc="best")
    a.set_xlim(-0.5, len(x) - 0.5)

save(fig, "figA_error_vs_frequency")

# ----------------------------------------------------------------------
# FIGURE C - spectrum, unmodulated against modulated
# ----------------------------------------------------------------------
print("Fig C")
if os.path.exists(FILE_UNMOD) and os.path.exists(FILE_PM):
    fig, ax = plt.subplots(2, 1, figsize=(7.2, 5.0), sharex=True)

    for a, path, title, col in (
            (ax[0], FILE_UNMOD, "Unmodulated carrier", C_HSE),
            (ax[1], FILE_PM, f"Phase modulated, "
                             r"$\beta$ = 1 rad", C_HSI)):
        t, v, dt = load(path)
        f, dB = spectrum(v, dt, ref_hz=F_CARRIER)
        a.plot(f / 1000, dB, color=col, lw=0.9)
        a.set_ylim(-80, 8)
        a.set_ylabel("level (dB rel. carrier)")
        a.set_title(title, loc="left")

    # mark the Bessel predictions on the modulated panel
    from math import factorial

    def bessel_j(n, x, terms=30):
        return sum((-1) ** k / (factorial(k) * factorial(k + n))
                   * (x / 2) ** (2 * k + n) for k in range(terms))

    j0 = bessel_j(0, BETA)
    for n in (1, 2, 3):
        lvl = 20 * np.log10(abs(bessel_j(n, BETA) / j0))
        for sgn in (-1, +1):
            fx = (F_CARRIER + sgn * n * F_MSG) / 1000
            ax[1].plot([fx], [lvl], marker="_", ms=14, mew=1.6,
                       color="k", zorder=5)
    ax[1].plot([], [], marker="_", ms=14, mew=1.6, color="k", ls="none",
               label=r"Bessel $J_n(\beta)/J_0(\beta)$")
    ax[1].legend(frameon=False, loc="upper right")

    ax[1].set_xlabel("frequency (kHz)")
    ax[1].set_xlim(2, 8)
    save(fig, "figC_spectrum_unmodulated_vs_pm")
else:
    print("  skipped: capture files not found")

# ----------------------------------------------------------------------
# FIGURE D - error by measurement point
# ----------------------------------------------------------------------
print("Fig D")
fig, ax = plt.subplots(figsize=(6.4, 3.2))
labels = [p[0] for p in POINTS]
means = np.array([p[1] for p in POINTS])
sds = np.array([p[2] for p in POINTS])
ns = [p[3] for p in POINTS]
x = np.arange(len(POINTS))

grand = means.mean()
ax.axhspan(grand - HSI_INSTABILITY, grand + HSI_INSTABILITY,
           color=C_REF, alpha=0.18,
           label=f"oscillator instability, $\\pm${HSI_INSTABILITY:.0f} ppm")
ax.axhline(grand, color=C_REF, ls="--", lw=1)
ax.errorbar(x, means, yerr=sds, fmt="o", color=C_HSI, capsize=4, ms=6,
            lw=1.3)
for i, (m, n) in enumerate(zip(means, ns)):
    ax.annotate(f"n = {n}", (i, m), textcoords="offset points",
                xytext=(9, 7), fontsize=7.5, color="0.35")

ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=8)
ax.set_ylabel("HSI frequency error (ppm)")
ax.set_title("Frequency error is independent of measurement point",
             loc="left")
ax.legend(frameon=False, loc="lower right")
ax.set_xlim(-0.5, len(x) - 0.5)
save(fig, "figD_error_by_measurement_point")

# ----------------------------------------------------------------------
# FIGURE E - cumulative timing error
# ----------------------------------------------------------------------
print("Fig E")
have = [p for p in (FILE_E_HSI, FILE_E_HSE) if os.path.exists(p)]
if have:
    fig, ax = plt.subplots(figsize=(6.4, 3.2))

    for path, label, col in ((FILE_E_HSI, "HSI", C_HSI),
                             (FILE_E_HSE, "HSE", C_HSE)):
        if not os.path.exists(path):
            continue
        t, v, dt = load(path)
        # rebuild an exact time axis, then find the zero crossings
        tt = t[0] + dt * np.arange(len(t))
        n = max(3, int(0.01 * len(v)))
        tt, vv = tt[n:-n], v[n:-n]
        vv = vv - vv.mean()
        thr = 0.25 * np.sqrt(2) * vv.std()
        raw = np.where((vv[:-1] < 0) & (vv[1:] >= 0))[0]
        keep, last = [], -1
        for i in raw:
            seg = vv[last + 1:i + 1]
            if seg.size and seg.min() < -thr:
                keep.append(i)
                last = i
        keep = np.array(keep)
        t0 = tt[keep] - vv[keep] * (tt[keep + 1] - tt[keep]) / \
            (vv[keep + 1] - vv[keep])
        k = np.round((t0 - t0[0]) / np.median(np.diff(t0)))
        ideal = t0[0] + k / F_CARRIER          # ideal times at exactly 5 kHz
        err = (t0 - ideal) * 1e6               # microseconds
        ax.plot((t0 - t0[0]) * 1000, err, color=col, lw=1.4, label=label)

    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel("elapsed time within the record (ms)")
    ax.set_ylabel("cumulative timing error (\u00b5s)")
    ax.set_title("Timing error accumulates at a rate set by the clock error",
                 loc="left")
    ax.legend(frameon=False)
    save(fig, "figE_cumulative_timing_error")
else:
    print("  skipped: capture files not found")

# ----------------------------------------------------------------------
# FIGURE F - uncertainty budget
# ----------------------------------------------------------------------
print("Fig F")
fig, ax = plt.subplots(1, 2, figsize=(7.6, 2.9), sharex=True)

for a, (clk, items) in zip(ax, BUDGET.items()):
    names = [n for n, _ in items][::-1]
    vals = [v for _, v in items][::-1]
    col = C_HSI if clk == "HSI" else C_HSE
    biggest = max(vals)
    colours = [col if v == biggest else "#bfc9ca" for v in vals]
    y = np.arange(len(vals))
    a.barh(y, vals, color=colours, height=0.6)
    for i, v in enumerate(vals):
        a.text(v * 1.15, i, f"{v:.1f}", va="center", fontsize=8)
    a.set_yticks(y)
    a.set_yticklabels(names, fontsize=8)
    a.set_xscale("log")
    a.set_xlim(0.5, 400)
    a.set_xlabel("contribution (ppm)")
    a.set_title(f"{clk}  -  limited by "
                + ("the oscillator" if clk == "HSI" else "the instrument"),
                loc="left", fontsize=9)
    a.grid(axis="y", visible=False)

save(fig, "figF_uncertainty_budget")

print(f"\nAll figures are in the '{OUT}' folder.")
print("Use the PDF versions in the thesis; they stay sharp at any size.")
