"""
m1_sweep_summary.py
=====================
Aggregates all M1 sweep captures into the headline result: fits the
realised phase sensitivity (slope), gain error, and linearity error
from commanded vs measured beta across the whole kp sweep.

Reuses extract_phase/measure_beta from check_m1_point.py.

Usage:
    Create a manifest text file, one capture per line: "kp  filename"

        0.25  m1_kp0p25_r1.csv
        0.25  m1_kp0p25_r2.csv
        0.50  m1_kp0p50_r1.csv
        ...

    python m1_sweep_summary.py m1_manifest.txt --fc 5000.12
"""

import argparse
import numpy as np
from check_m1_point import extract_phase, measure_beta, TRIM_CYCLES, BW_HZ, FM_HZ
from measure_fc_ref import load_tbs_csv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("manifest")
    ap.add_argument("--fc", type=float, required=True)
    ap.add_argument("--bw", type=float, default=BW_HZ)
    args = ap.parse_args()

    rows = []
    with open(args.manifest) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            kp_str, fname = line.split(None, 1)
            kp = float(kp_str)
            fname = fname.strip()
            t, x = load_tbs_csv(fname, channel=1)
            t_p, theta, _ = extract_phase(t, x, args.fc, args.bw, TRIM_CYCLES)
            beta_hat = measure_beta(t_p, theta, FM_HZ)
            rows.append((kp, beta_hat, fname))
            ge = (beta_hat / kp - 1) * 100 if kp != 0 else float("nan")
            print(f"kp={kp:.2f}  beta_hat={beta_hat:.4f}  "
                  f"gain_err={ge:+.2f}%  [{fname}]")

    kps = np.array([r[0] for r in rows])
    betas = np.array([r[1] for r in rows])

    A = np.vstack([kps, np.ones_like(kps)]).T
    coef, *_ = np.linalg.lstsq(A, betas, rcond=None)
    slope, intercept = coef
    pred = A @ coef
    ss_res = np.sum((betas - pred) ** 2)
    ss_tot = np.sum((betas - betas.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot

    resid = betas - pred
    linearity_pct = float(np.max(np.abs(resid)) / 4.0 * 100)  # % of 4.0 rad full scale
    gain_error_pct = float((slope - 1) * 100)

    print(f"\n=== M1 sweep fit ({len(rows)} points) ===")
    print(f"Slope (realised phase sensitivity) : {slope:.4f}")
    print(f"Gain error                          : {gain_error_pct:+.2f} %")
    print(f"Intercept                           : {intercept:+.4f} rad")
    print(f"R^2                                  : {r2:.6f}")
    print(f"Max linearity deviation              : {linearity_pct:.2f} % of full scale")

    if r2 < 0.999:
        print("\nNOTE: R^2 below 0.999 -- check for an outlier point before "
              "trusting the slope.")


if __name__ == "__main__":
    main()
