"""
plot_drift.py   -   T5 warm-up drift

Reads every T5 capture in a folder, works out the elapsed time from each
file name, and plots frequency drift against time for HSI and HSE.

FILE NAMING
-----------
Put the elapsed time in seconds into the name, like this:

    T5_HSI_t0000s.csv     0 seconds
    T5_HSI_t0030s.csv     30 seconds
    T5_HSI_t0300s.csv     5 minutes
    T5_HSE_t1800s.csv     30 minutes

Minutes also work: T5_HSI_t05min.csv

WHY THIS RESULT IS NOT INSTRUMENT-LIMITED
-----------------------------------------
Drift is measured against the instrument's own earlier reading, so the
oscilloscope timebase error is common to every point and cancels. Unlike
absolute accuracy, this measurement is limited only by the estimator,
which is about 1 ppm here.
"""

import os
import re
import hashlib
import numpy as np

from analyse_freq import load_tek_csv, freq_zero_crossing, ppm

# ---- settings ---------------------------------------------------------
FOLDER = "."
F_SET = 5000.0
SETTLE_BAND_PPM = 5.0      # drift is "settled" once inside this band
N_REF = 3                  # how many final captures define the settled value
# -----------------------------------------------------------------------


def elapsed_seconds(name):
    """Read the elapsed time out of a file name."""
    m = re.search(r"t?(\d+(?:\.\d+)?)\s*min", name, re.I)
    if m:
        return float(m.group(1)) * 60.0
    m = re.search(r"t(\d+(?:\.\d+)?)\s*s", name, re.I)
    if m:
        return float(m.group(1))
    m = re.search(r"(\d+(?:\.\d+)?)\s*s\b", name, re.I)
    if m:
        return float(m.group(1))
    return None


def clock_of(name):
    up = name.upper()
    return "HSI" if "HSI" in up else ("HSE" if "HSE" in up else "?")


def file_hash(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


files = sorted(f for f in os.listdir(FOLDER)
               if f.lower().endswith(".csv") and "result" not in f.lower())
if not files:
    print(f"No CSV files in {os.path.abspath(FOLDER)}")
    raise SystemExit(1)

print(f"Folder: {os.path.abspath(FOLDER)}")
print(f"{len(files)} CSV file(s)\n")

# ---- duplicate check --------------------------------------------------
# Saving twice without pressing Single first returns the same frozen record.
# A duplicate looks like a perfectly flat section of the drift curve, so it
# has to be caught before plotting.
seen = {}
dupes = []
for name in files:
    h = file_hash(os.path.join(FOLDER, name))
    if h in seen:
        dupes.append((name, seen[h]))
    else:
        seen[h] = name

if dupes:
    print("DUPLICATE FILES FOUND - these are byte-identical:")
    for a, b in dupes:
        print(f"   {a}  is a copy of  {b}")
    print("   Cause: saving without pressing Single first, so the scope")
    print("   wrote out the same frozen record twice. Exclude or recapture")
    print("   these before reporting a drift curve.\n")
else:
    print("No duplicate files. Every capture is distinct.\n")

# ---- analyse ----------------------------------------------------------
data = {"HSI": [], "HSE": []}
skipped = []

for name in files:
    t_s = elapsed_seconds(name)
    clk = clock_of(name)
    if t_s is None or clk == "?":
        skipped.append((name, "cannot read time or clock from the name"))
        continue
    try:
        d = load_tek_csv(os.path.join(FOLDER, name))
        r = freq_zero_crossing(d["t"], d["channels"][0], dt=d["dt"])
        if not r["ok"]:
            skipped.append((name, f"status CHECK: {r['notes']}"))
            continue
        data[clk].append((t_s, r["freq"], name))
    except Exception as e:
        skipped.append((name, str(e)))

for clk in data:
    data[clk].sort(key=lambda x: x[0])

if skipped:
    print(f"{len(skipped)} file(s) skipped:")
    for n, why in skipped:
        print(f"   {n}: {why}")
    print()

# ---- drift relative to the settled value ------------------------------
results = {}
for clk, rows in data.items():
    if len(rows) < 3:
        continue
    t = np.array([r[0] for r in rows])
    f = np.array([r[1] for r in rows])
    f_ref = f[-N_REF:].mean()
    drift = (f - f_ref) / f_ref * 1e6
    results[clk] = dict(t=t, f=f, drift=drift, f_ref=f_ref,
                        names=[r[2] for r in rows])

    print(f"--- {clk} ---")
    print(f"  captures        : {len(t)} over {t[-1]/60:.1f} minutes")
    print(f"  settled value   : {f_ref:.4f} Hz "
          f"({ppm(f_ref, F_SET):+.1f} ppm from target)")
    print(f"  drift at t = 0  : {drift[0]:+.2f} ppm")
    print(f"  total excursion : {drift.max()-drift.min():.2f} ppm")

    inside = np.abs(drift) < SETTLE_BAND_PPM
    settle = None
    for i in range(len(t)):
        if inside[i:].all():
            settle = t[i]
            break
    if settle is None:
        print(f"  settling time   : not settled within "
              f"{SETTLE_BAND_PPM:.0f} ppm during the run")
    elif settle == 0:
        print(f"  settling time   : already inside {SETTLE_BAND_PPM:.0f} ppm "
              f"at t = 0 (no measurable warm-up drift)")
    else:
        print(f"  settling time   : {settle:.0f} s "
              f"({settle/60:.1f} min) to stay inside "
              f"{SETTLE_BAND_PPM:.0f} ppm")
    print()

# ---- write the table --------------------------------------------------
out = os.path.join(FOLDER, "T5_drift_results.csv")
with open(out, "w") as fh:
    fh.write("clock,elapsed_s,elapsed_min,f_measured_Hz,"
             "drift_vs_settled_ppm,error_vs_target_ppm,file\n")
    for clk, r in results.items():
        for i in range(len(r["t"])):
            fh.write(f"{clk},{r['t'][i]:.0f},{r['t'][i]/60:.2f},"
                     f"{r['f'][i]:.4f},{r['drift'][i]:.3f},"
                     f"{ppm(r['f'][i], F_SET):.2f},{r['names'][i]}\n")
print(f"Table written to {out}")

# ---- plot -------------------------------------------------------------
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    style = {"HSI": dict(color="#c0392b", marker="o"),
             "HSE": dict(color="#2874a6", marker="s")}

    for clk, r in results.items():
        st = style.get(clk, {})
        ax[0].plot(r["t"] / 60, r["drift"], ms=4, lw=1.2, label=clk, **st)
        ax[1].plot(r["t"] / 60, r["drift"], ms=4, lw=1.2, label=clk, **st)

    ax[0].set_title("Warm-up drift, both clock sources")
    ax[1].set_title("Same data, zoomed to the HSE scale")
    if "HSE" in results:
        m = np.abs(results["HSE"]["drift"]).max()
        ax[1].set_ylim(-max(m * 1.6, 3), max(m * 1.6, 3))
    for a in ax:
        a.axhline(0, color="k", lw=0.6)
        a.axhspan(-SETTLE_BAND_PPM, SETTLE_BAND_PPM, color="grey", alpha=0.12)
        a.set_xlabel("elapsed time (minutes)")
        a.set_ylabel("drift from settled value (ppm)")
        a.grid(alpha=0.3)
        a.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(FOLDER, "T5_warmup_drift.png"), dpi=150)
    fig.savefig(os.path.join(FOLDER, "T5_warmup_drift.pdf"))
    print("Plots written to T5_warmup_drift.png and .pdf")
except Exception as e:
    print(f"(no plot: {e})")
