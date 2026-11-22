"""
m3_simple_demo.py
==================
Produces the "simple oscilloscope version" figure for Section 6.4:
the BPSK phase reversal, shown as a time shift between zero crossings.

No separate unmodulated-carrier capture is needed. Symbol 0 of the bit
pattern (bits[0] = 0) carries no added phase, so it IS the unmodulated
carrier -- generated through the same code path, same sample rate, same
amplitude as every other symbol. This script compares the zero-crossing
timing of that first ("0") segment against a later ("1") segment of the
SAME capture.

Usage:
    python m3_simple_demo.py m3_bpsk_r1.csv --fc 5000.123

    --fc : measured carrier frequency from R1 (fc_ref.json). Only used to
           group crossings into cycles -- the result does not depend on
           its precision to more than ~100 ppm.

Cross-check: the Delta-phi printed here should agree with the "Separation"
column of Table 10.4 (M3 quantitative Hilbert analysis) to within a few
degrees. If they disagree by more than that, something is wrong with one
of the two analyses -- investigate before trusting either.
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt

RS = 250.0                 # BPSK symbol rate, sym/s
T_SYM = 1.0 / RS            # 4 ms
GUARD = 400e-6              # exclude this much around every transition
bits = [0, 1, 0, 1, 1, 0, 1, 0, 0, 1]


# ----------------------------------------------------------------------
def load_tbs_csv(path, channel=1):
    """Robust TBS1000C CSV loader (same as measure_fc_ref.py)."""
    rows = []
    with open(path, "r", errors="ignore") as f:
        for line in f:
            parts = [p.strip() for p in line.split(",")]
            vals = []
            for p in parts:
                try:
                    vals.append(float(p))
                except ValueError:
                    vals.append(np.nan)
            if np.count_nonzero(~np.isnan(vals)) >= 2:
                rows.append(vals)
    width = max(len(r) for r in rows)
    data = np.full((len(rows), width), np.nan)
    for i, r in enumerate(rows):
        data[i, :len(r)] = r
    good = np.mean(~np.isnan(data), axis=0) > 0.95
    data = data[:, good]
    data = data[~np.isnan(data).any(axis=1)]

    tcol = None
    for c in range(data.shape[1]):
        d = np.diff(data[:, c])
        if np.all(d > 0) and np.std(d) / np.mean(d) < 1e-3:
            tcol = c
            break
    if tcol is None:
        raise ValueError("Could not identify the time column")
    vcols = [c for c in range(data.shape[1]) if c != tcol]
    return data[:, tcol], data[:, vcols[channel - 1]]


def zero_crossings(t, x, hyst_frac=0.15):
    """Rising zero crossings with Schmitt hysteresis, sub-sample interpolated."""
    xa = x - np.median(x)
    amp = (np.percentile(xa, 99.5) - np.percentile(xa, 0.5)) / 2.0
    h = hyst_frac * amp
    t0, armed = [], False
    for i in range(len(xa) - 1):
        if xa[i] < -h:
            armed = True
        if armed and xa[i] < 0.0 <= xa[i + 1]:
            frac = -xa[i] / (xa[i + 1] - xa[i])
            t0.append(t[i] + frac * (t[i + 1] - t[i]))
            armed = False
    return np.asarray(t0)


# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--channel", type=int, default=1)
    ap.add_argument("--fc", type=float, required=True,
                    help="measured carrier freq from R1, e.g. 5000.12")
    args = ap.parse_args()

    t, x = load_tbs_csv(args.csv, args.channel)
    zc = zero_crossings(t, x)
    Tc = 1.0 / args.fc

    # ---- reference trend: symbol 0, which carries bits[0] = 0 (unmodulated) ----
    mask_ref = (zc > GUARD) & (zc < T_SYM - GUARD)
    zc_ref = zc[mask_ref]
    if len(zc_ref) < 4:
        raise ValueError("too few crossings in the reference segment -- "
                          "check --fc and the guard interval")

    k_ref = np.round((zc_ref - zc_ref[0]) / Tc)
    A = np.vstack([k_ref, np.ones_like(k_ref)]).T
    coef, *_ = np.linalg.lstsq(A, zc_ref, rcond=None)  # [period, t0]
    predict = lambda k: coef[0] * k + coef[1]

    # ---- find the first bit-transition boundary and pick the modulated
    #      segment right after it ----
    trans_idx = next(i for i in range(1, len(bits)) if bits[i] != bits[i - 1])
    t_trans = trans_idx * T_SYM
    mask_mod = (zc > t_trans + GUARD) & (zc < t_trans + T_SYM - GUARD)
    zc_mod = zc[mask_mod]
    if len(zc_mod) < 4:
        raise ValueError("too few crossings in the modulated segment")

    k_mod = np.round((zc_mod - coef[1]) / coef[0])
    predicted = predict(k_mod)
    delta_t = float(np.mean(zc_mod - predicted))
    delta_phi = 360.0 * delta_t / Tc
    # wrap to the +/-180 range for a sane printed value
    delta_phi_wrapped = ((delta_phi + 180.0) % 360.0) - 180.0

    print(f"Reference segment : bit index 0 (bits[0] = {bits[0]}), "
          f"{len(zc_ref)} crossings")
    print(f"Modulated segment : bit index {trans_idx} "
          f"(bits[{trans_idx}] = {bits[trans_idx]}), "
          f"transition at t = {t_trans * 1e3:.3f} ms, "
          f"{len(zc_mod)} crossings")
    print(f"Delta t   = {delta_t * 1e6:+.2f} us")
    print(f"Delta phi = {delta_phi_wrapped:+.1f} deg   (ideal: 180.0 deg)")
    print("Compare this against the 'Separation' column of Table 10.4 "
          "(M3 Hilbert analysis) -- they should agree to within a few "
          "degrees.")

    # ---------------------------------------------------------------
    # Figure: zoom into the transition, show the ideal (dashed) trace
    # continuing straight through against the actual signal bending away.
    # ---------------------------------------------------------------
    span = 5 * Tc  # a few carrier cycles either side of the transition
    lo, hi = t_trans - span, t_trans + span
    sl = (t >= lo) & (t <= hi)

    tt = np.linspace(lo, hi, 2000)
    ideal = np.interp(tt, zc_ref, np.zeros_like(zc_ref))  # placeholder, see below
    # build a clean ideal cosine using the fitted period/phase instead:
    phase0 = -2 * np.pi * coef[1] / coef[0]  # phase at t=0 implied by the fit
    ideal_wave = np.cos(2 * np.pi * tt / coef[0] + phase0)
    # scale ideal to the measured amplitude for a fair overlay
    amp = (np.percentile(x - np.median(x), 99.5)
           - np.percentile(x - np.median(x), 0.5)) / 2.0
    ideal_wave = ideal_wave * amp + np.median(x)

    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(t[sl] * 1e3, x[sl], color="tab:blue", lw=1.5, label="Measured (LPF output)")
    ax.plot(tt * 1e3, ideal_wave, "--", color="tab:gray", lw=1.2,
            label="Unmodulated carrier (bit 0 trend, extrapolated)")
    ax.axvline(t_trans * 1e3, color="tab:red", ls=":", lw=1,
               label=f"Commanded transition ({bits[trans_idx-1]}\u2192{bits[trans_idx]})")

    for zc_t in zc_ref[zc_ref > t_trans - span]:
        ax.axvline(zc_t * 1e3, color="tab:gray", alpha=0.25, lw=0.8)
    for zc_t in zc_mod[zc_mod < t_trans + span]:
        ax.axvline(zc_t * 1e3, color="tab:blue", alpha=0.35, lw=0.8)

    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("LPF output (V)")
    ax.set_title(f"BPSK phase reversal: measured $\\Delta t$ = "
                 f"{delta_t*1e6:+.1f} $\\mu$s "
                 f"($\\Delta\\phi$ = {delta_phi_wrapped:+.1f}$^\\circ$)")
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig("m3_phase_reversal_demo.png", dpi=200)
    print("\nFigure saved: m3_phase_reversal_demo.png")


if __name__ == "__main__":
    main()
