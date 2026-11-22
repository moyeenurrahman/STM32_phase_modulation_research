"""
batch_freq.py

Analyses every CSV in a folder in one run and writes a results table.

The target frequency is read from the file name, so name your files with
the frequency in them, for example:

    T1_HSI_1000Hz_run1.csv
    T1_HSI_10kHz_run3.csv
    T4_HSI_8MHz_run2.csv
    T6_HSE_5000Hz_run1.csv

Anything containing "T6" or "PM" is analysed with Method C (interpolated
FFT) because modulation displaces the zero crossings. Everything else uses
Method D (zero-crossing regression).

Press Run in PyCharm. Results are printed and also written to
frequency_results.csv, which opens directly in Excel.
"""

import os
import re
import numpy as np

from analyse_freq import (load_tek_csv, freq_zero_crossing, freq_fft_interp,
                          is_sample_and_hold, ppm)

# ---- settings ---------------------------------------------------------
FOLDER = "."              # "." means the project folder
DEFAULT_F = 5000.0        # used if the file name has no frequency in it
SCOPE_TIMEBASE_PPM = 20.0
# -----------------------------------------------------------------------


def freq_from_name(name):
    """Pull the target frequency out of a file name."""
    m = re.search(r"(\d+(?:\.\d+)?)\s*([kKmM]?)\s*[hH][zZ]", name)
    if not m:
        return None
    value = float(m.group(1))
    prefix = m.group(2).lower()
    if prefix == "k":
        value *= 1e3
    elif prefix == "m":
        value *= 1e6
    return value


def test_and_clock(name):
    """Pull the test id (T1..T7) and clock source out of the file name."""
    up = name.upper()
    m = re.search(r"\bT([1-7])\b|T([1-7])[_\-]", up)
    test = f"T{m.group(1) or m.group(2)}" if m else "?"
    if "HSI" in up:
        clk = "HSI"
    elif "HSE" in up:
        clk = "HSE"
    else:
        clk = "?"
    return test, clk


def is_pm_file(name):
    upper = name.upper()
    return "T6" in upper or "_PM" in upper or "PM_" in upper


def analyse_one(path, f_set, pm, ch=0, d=None):
    if d is None:
        d = load_tek_csv(path)
    v = d["channels"][ch]
    dt = d["dt"]

    if pm:
        r = freq_fft_interp(d["t"], v, dt=dt)
        method = "C"
        sah = False
    else:
        r = freq_zero_crossing(d["t"], v, dt=dt)
        method = "D"
        sah = is_sample_and_hold(v, dt, f_set)

    return dict(freq=r["freq"], err=ppm(r["freq"], f_set), method=method,
                sigma=r.get("sigma_ppm", float("nan")),
                jitter_ns=r.get("jitter_rms", 0.0) * 1e9,
                ok=r["ok"], sah=sah, notes=r["notes"],
                vpp=v.max() - v.min(), dt=dt, npts=len(v))


files = sorted(f for f in os.listdir(FOLDER)
               if f.lower().endswith(".csv") and f != "frequency_results.csv")

if not files:
    print(f"No CSV files found in {os.path.abspath(FOLDER)}")
    raise SystemExit(1)

print(f"Folder: {os.path.abspath(FOLDER)}")
print(f"{len(files)} CSV file(s) found\n")

rows = []
for name in files:
    path = os.path.join(FOLDER, name)
    f_set = freq_from_name(name)
    if f_set is None:
        f_set = DEFAULT_F
        guessed = " (frequency guessed)"
    else:
        guessed = ""
    pm = is_pm_file(name)

    try:
        d = load_tek_csv(path)
        nch = len(d["channels"])
        for c in range(nch):
            r = analyse_one(path, f_set, pm, ch=c, d=d)
            label = name if nch == 1 else f"{name} [ch{c + 1}]"
            r["file"] = label
            r["f_set"] = f_set
            r["test"], r["clock"] = test_and_clock(name)
            rows.append(r)
            flag = ""
            if not r["ok"]:
                flag = "  <-- CHECK"
            elif r["sah"]:
                flag = "  <-- staircase, sigma not valid"
            print(f"{label:38s} {f_set:10.0f} Hz  {r['err']:+10.2f} ppm  "
                  f"[{r['method']}]{flag}{guessed}")
    except Exception as e:
        print(f"{name:38s} FAILED: {e}")
        tst, clk = test_and_clock(name)
        rows.append(dict(file=name, f_set=f_set, test=tst, clock=clk,
                         freq=float("nan"),
                         err=float("nan"), method="-", sigma=float("nan"),
                         jitter_ns=float("nan"), ok=False, sah=False,
                         notes=str(e), vpp=float("nan"), dt=float("nan"),
                         npts=0))

# ---- write the results table -----------------------------------------
out = os.path.join(FOLDER, "frequency_results.csv")
with open(out, "w") as f:
    f.write("file,test,clock,f_set_Hz,f_measured_Hz,abs_error_Hz,"
            "signed_error_ppm,method,estimator_sigma_ppm,jitter_ns,Vpp,"
            "points,status,notes\n")
    for r in rows:
        status = "OK" if r["ok"] else "CHECK"
        if r["ok"] and r["sah"]:
            status = "STAIRCASE"
        f.write(f"{r['file']},{r.get('test','?')},{r.get('clock','?')},"
                f"{r['f_set']:.1f},{r['freq']:.4f},"
                f"{r['freq'] - r['f_set']:.4f},{r['err']:.2f},{r['method']},"
                f"{r['sigma']:.3f},{r['jitter_ns']:.0f},{r['vpp']:.3f},"
                f"{r['npts']},{status},{r['notes']}\n")
print(f"\nWritten to {out}")

# ---- group by target frequency ---------------------------------------
good = [r for r in rows if r["ok"] and not r["sah"] and np.isfinite(r["err"])]
if good:
    print("\n--- Grouped by test, clock and frequency ---")
    print(f"{'test':>5} {'clock':>6} {'f_set':>10}  {'runs':>4}  "
          f"{'mean ppm':>11}  {'1 sigma':>9}")
    keys = sorted({(r["test"], r["clock"], r["f_set"]) for r in good})
    for k in keys:
        e = np.array([r["err"] for r in good
                      if (r["test"], r["clock"], r["f_set"]) == k])
        sd = e.std(ddof=1) if len(e) > 1 else 0.0
        print(f"{k[0]:>5} {k[1]:>6} {k[2]:10.0f}  {len(e):4d}  "
              f"{e.mean():+11.2f}  {sd:9.2f}")

    # Consistency across frequencies, done separately for each clock source.
    # All frequencies share one oscillator, so they must agree in sign and
    # magnitude. A large spread means a firmware fault, not a clock error.
    for clk in sorted({r["clock"] for r in good if r["clock"] != "?"}):
        sub = [r for r in good if r["clock"] == clk and r["test"] in
               ("T1", "T2")]
        means = {fs: np.mean([r["err"] for r in sub if r["f_set"] == fs])
                 for fs in sorted({r["f_set"] for r in sub})}
        if len(means) < 2:
            continue
        spread = max(means.values()) - min(means.values())
        print(f"\n{clk} consistency across {len(means)} frequencies: "
              f"spread {spread:.2f} ppm")
        if len({np.sign(v) for v in means.values()}) > 1:
            print("  >> The frequencies do not share a common sign. Check the")
            print("     firmware cycle counts before reporting.")
        elif spread > 10 * SCOPE_TIMEBASE_PPM:
            print("  >> Spread is large. Check PSC, ARR and the cycle count")
            print("     in the DMA buffer for each frequency.")
        else:
            print("  All frequencies agree in sign and magnitude, which is")
            print("  what an error originating in the clock source looks")
            print("  like. The generation chain is consistent.")

    # T4 versus the DAC path, the strongest result available
    for clk in ("HSI", "HSE"):
        mco = [r["err"] for r in good if r["test"] == "T4" and
               r["clock"] == clk]
        dac = [r["err"] for r in good if r["test"] in ("T1", "T2") and
               r["clock"] == clk]
        if mco and dac:
            print(f"\n{clk}: MCO pin {np.mean(mco):+.2f} ppm vs DAC path "
                  f"{np.mean(dac):+.2f} ppm  "
                  f"(difference {np.mean(mco) - np.mean(dac):+.2f} ppm)")
            print("  Agreement here shows the DAC, DMA and timer add no")
            print("  frequency error beyond that of the oscillator.")

bad = [r for r in rows if not r["ok"]]
if bad:
    print(f"\n{len(bad)} file(s) need attention:")
    for r in bad:
        print(f"   {r['file']}: {r['notes']}")
