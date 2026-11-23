"""
diagnose_csv_load.py
======================
Standalone diagnostic for a CSV that load_tbs_csv fails to parse.
Prints the raw header lines and, for every numeric column found, the
statistics used to decide whether it's the time column -- so we can
see exactly why detection is failing, instead of guessing further.

Usage:
    python diagnose_csv_load.py m3_transition.csv
"""

import sys
import numpy as np

path = sys.argv[1]

print("=== First 15 raw lines ===")
raw_lines = []
with open(path, "r", errors="ignore") as f:
    for i, line in enumerate(f):
        raw_lines.append(line.rstrip("\n"))
        if i < 15:
            print(f"{i:3d}: {line.rstrip()}")

print(f"\nTotal lines in file: {len(raw_lines)}")

rows = []
for line in raw_lines:
    parts = [p.strip() for p in line.split(",")]
    vals = []
    for p in parts:
        try:
            vals.append(float(p))
        except ValueError:
            vals.append(np.nan)
    if np.count_nonzero(~np.isnan(vals)) >= 2:
        rows.append(vals)

print(f"Rows kept as 'data rows' (>=2 numeric fields): {len(rows)}")

if not rows:
    print("\nNo rows had >=2 numeric fields -- the file may use a "
          "different delimiter, or every field failed to parse as a "
          "plain float. Check the raw lines above for anything unusual "
          "(semicolons, units attached to numbers, stray quotes).")
    sys.exit(0)

width = max(len(r) for r in rows)
data = np.full((len(rows), width), np.nan)
for i, r in enumerate(rows):
    data[i, :len(r)] = r

good = np.mean(~np.isnan(data), axis=0) > 0.95
data_g = data[:, good]
print(f"Columns before/after >95%-numeric filter: "
      f"{data.shape[1]} -> {data_g.shape[1]}")

data_g = data_g[~np.isnan(data_g).any(axis=1)]
print(f"Rows after dropping any remaining NaNs: {data_g.shape[0]}")

print(f"\n=== Per-column diagnostics ({data_g.shape[1]} columns) ===")
for c in range(data_g.shape[1]):
    col = data_g[:, c]
    d = np.diff(col)
    if len(d) == 0:
        print(f"col {c}: no data")
        continue
    frac_positive = float(np.mean(d > 0))
    med = float(np.median(d))
    frac_near_median = (float(np.mean(np.abs(d - med) < 0.05 * abs(med)))
                        if med != 0 else float("nan"))
    print(f"col {c}: first={col[0]:.6g}  last={col[-1]:.6g}  "
          f"median_step={med:.6g}  frac_positive={frac_positive:.3f}  "
          f"frac_near_median={frac_near_median:.3f}")
