"""
analyse_stability.py   -   T5 reframed as short-term frequency stability

Reads T5_drift_results.csv (written by plot_drift.py) and asks the questions
the data can actually answer:

  1. Is there a statistically significant drift trend, or only scatter?
  2. How stable is each oscillator over minutes?
  3. How much of the scatter is the measurement system rather than the
     oscillator?
  4. If no trend is resolved, what upper bound can be placed on drift?

Put this next to T5_drift_results.csv and press Run.
"""

import csv
import numpy as np
from math import lgamma

FILE = "T5_drift_results.csv"
ESTIMATOR_PPM = 1.1        # from the Session 0 validation


def t_pvalue(t, df):
    """Two-sided p-value for Student's t, by numerical integration."""
    x = df / (df + t * t)
    a, b = df / 2.0, 0.5
    u = np.linspace(1e-12, x, 400000)
    I = np.trapezoid(u ** (a - 1) * (1 - u) ** (b - 1), u)
    I *= np.exp(lgamma(a + b) - lgamma(a) - lgamma(b))
    return float(np.clip(I, 0.0, 1.0))


rows = list(csv.DictReader(open(FILE)))
data = {}
for r in rows:
    data.setdefault(r["clock"], []).append(
        (float(r["elapsed_min"]), float(r["drift_vs_settled_ppm"])))

summary = {}
for clk in sorted(data):
    a = np.array(sorted(data[clk]))
    t, y = a[:, 0], a[:, 1]
    n = len(t)

    # ---- linear trend test -------------------------------------------
    slope, intercept = np.polyfit(t, y, 1)
    resid = y - (slope * t + intercept)
    s_resid = resid.std(ddof=2)
    se_slope = s_resid / np.sqrt(((t - t.mean()) ** 2).sum())
    tstat = slope / se_slope
    p = t_pvalue(tstat, n - 2)

    # ---- scatter and stability ---------------------------------------
    sd = y.std(ddof=1)

    # ---- Allan deviation at the sampling interval ---------------------
    # Differences between consecutive readings measure short-term
    # stability without being affected by any slow trend.
    d = np.diff(y)
    allan = np.sqrt(np.mean(d ** 2) / 2.0)

    summary[clk] = dict(n=n, sd=sd, slope=slope, p=p, allan=allan,
                        span=t[-1], rng=y.max() - y.min(),
                        se_slope=se_slope)

    print(f"--- {clk} ---")
    print(f"  captures                 : {n} over {t[-1]:.1f} min")
    print(f"  scatter, 1 sigma         : {sd:8.2f} ppm")
    print(f"  peak to peak             : {y.max()-y.min():8.2f} ppm")
    print(f"  short-term stability     : {allan:8.2f} ppm  "
          f"(Allan dev, consecutive pairs)")
    print(f"  trend                    : {slope:+8.3f} ppm/min "
          f"+/- {se_slope:.3f}")
    print(f"  trend significance       : t = {tstat:+.2f}, p = {p:.3f}  "
          f"-> {'SIGNIFICANT' if p < 0.05 else 'not significant'}")
    if p >= 0.05:
        bound = 2 * se_slope * t[-1]
        print(f"  upper bound on total drift: |drift| < {bound:.0f} ppm "
              f"over {t[-1]:.0f} min (95 % confidence)")
    print()

# ---- compare the two clocks ------------------------------------------
if "HSI" in summary and "HSE" in summary:
    hi, he = summary["HSI"], summary["HSE"]
    F = (hi["sd"] ** 2) / (he["sd"] ** 2)
    print("--- HSI compared with HSE ---")
    print(f"  scatter ratio            : {hi['sd']/he['sd']:.1f}x "
          f"(variance ratio F = {F:.1f})")
    print(f"  HSI is {hi['sd']/he['sd']:.0f} times less stable than the "
          f"crystal reference.")
    print()

    # HSE scatter sets the floor for the whole measurement: a crystal
    # cannot really wander by tens of ppm in minutes, so that scatter is
    # instrumental. Subtracting it in quadrature isolates the HSI part.
    print("--- separating oscillator from instrument ---")
    print(f"  HSE scatter              : {he['sd']:8.2f} ppm")
    print(f"    estimator contribution : {ESTIMATOR_PPM:8.2f} ppm")
    inst = np.sqrt(max(he["sd"] ** 2 - ESTIMATOR_PPM ** 2, 0))
    print(f"    remaining, capture to")
    print(f"    capture instrumental   : {inst:8.2f} ppm")
    print("    A quartz crystal cannot wander this much in minutes, so")
    print("    this is taken as the measurement floor.")
    print()
    intrinsic = np.sqrt(max(hi["sd"] ** 2 - he["sd"] ** 2, 0))
    print(f"  HSI scatter              : {hi['sd']:8.2f} ppm")
    print(f"    minus the floor above  : {intrinsic:8.2f} ppm")
    print("    This is the intrinsic short-term instability of the RC")
    print("    oscillator, over and above the measurement floor.")
