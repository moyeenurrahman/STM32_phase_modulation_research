"""
m3_transition_detail.py
=========================
Analyzes the M3 zoomed transition capture: locates the phase transition
automatically, computes the 10-90% transition time and the peak
overshoot beyond the settled level -- the last two columns of Table
10.4 that the whole-buffer analysis (m3_sweep_summary.py) doesn't cover.

Uses a WIDER post-detection bandwidth than the steady-state RMS
analysis (default 15 kHz vs 2.5 kHz). A narrow-band filter here would
artificially slow the observed transition with its own settling on top
of the real analog filter's dynamics, biasing the transition-time
measurement upward. 15 kHz sits just under the LPF's ~20.1 kHz cutoff,
so it captures the true step response shape without adding extra
smoothing.

Usage:
    python m3_transition_detail.py m3_transition.csv --fc 5000.12

Produces m3_transition_detail.png marking the detected transition and
the 10%/90% crossing points. CHECK THIS PLOT before trusting the
printed numbers -- automatic transition detection can occasionally
lock onto the wrong feature, especially if the capture wasn't well
centered on the transition.
"""

import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import hilbert, butter, filtfilt
from measure_fc_ref import load_tbs_csv

BW_HZ = 15000.0     # wide -- NOT the same 2500 Hz used for steady-state RMS
TRIM_CYCLES = 3      # fewer cycles trimmed -- this record is much shorter
                     # than the wide M1/M3 captures


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--fc", type=float, required=True)
    ap.add_argument("--bw", type=float, default=BW_HZ)
    ap.add_argument("--trim-cycles", type=int, default=TRIM_CYCLES)
    ap.add_argument("--settle-window-us", type=float, default=300.0,
                    help="window on each side used to establish the "
                         "before/after phase levels")
    ap.add_argument("--overshoot-window-us", type=float, default=400.0,
                    help="how far past the transition to search for "
                         "the overshoot peak")
    args = ap.parse_args()

    t, x = load_tbs_csv(args.csv, channel=1)
    fs = 1.0 / ((t[-1] - t[0]) / (len(t) - 1))
    dt = 1.0 / fs
    print(f"record length  : {len(t)}")
    print(f"sample interval : {dt*1e6:.4f} us")
    print(f"record span     : {(t[-1]-t[0])*1e3:.3f} ms")

    x_ac = x - np.mean(x)
    z = hilbert(x_ac)
    zb = z * np.exp(-1j * 2 * np.pi * args.fc * t)
    b, a = butter(4, args.bw / (fs / 2.0))
    zb = filtfilt(b, a, zb.real) + 1j * filtfilt(b, a, zb.imag)
    theta = np.unwrap(np.angle(zb))
    env = np.abs(zb)

    n = int(args.trim_cycles * fs / args.fc)
    if 2 * n >= len(t):
        raise ValueError("trim removes the entire record -- reduce "
                          "--trim-cycles or check --fc")
    sl = slice(n, len(t) - n)
    t_p, theta_p, env_p = t[sl], theta[sl], env[sl]

    # locate the transition: point of steepest change in unwrapped phase
    dtheta = np.abs(np.diff(theta_p))
    i_trans = int(np.argmax(dtheta))
    print(f"\nApproximate transition at t = {t_p[i_trans]*1e3:.4f} ms")

    settle_n = int(args.settle_window_us * 1e-6 / dt)
    if i_trans - settle_n < 0 or i_trans + settle_n >= len(t_p):
        print("WARNING: not enough margin on one side to establish the "
              "settled levels -- widen the capture or reduce "
              "--settle-window-us")

    before = theta_p[max(0, i_trans - settle_n):i_trans]
    guard_n = int(50e-6 / dt)   # skip a short guard after the jump before
                                # sampling the "after" level
    after = theta_p[i_trans + guard_n: i_trans + guard_n + settle_n]

    before_level = float(np.median(before))
    after_level = float(np.median(after))
    step = after_level - before_level

    print(f"Before level : {np.degrees(before_level):+.2f} deg")
    print(f"After level  : {np.degrees(after_level):+.2f} deg")
    print(f"Phase step   : {np.degrees(step):+.2f} deg")

    level_10 = before_level + 0.1 * step
    level_90 = before_level + 0.9 * step

    search_n = int(500e-6 / dt)
    lo = max(0, i_trans - search_n)
    hi = min(len(t_p), i_trans + search_n)
    seg_t = t_p[lo:hi]
    seg_th = theta_p[lo:hi]

    def first_crossing(level, rising):
        if rising:
            idx = np.where((seg_th[:-1] < level) & (seg_th[1:] >= level))[0]
        else:
            idx = np.where((seg_th[:-1] > level) & (seg_th[1:] <= level))[0]
        if len(idx) == 0:
            return None
        i = int(idx[0])
        frac = (level - seg_th[i]) / (seg_th[i + 1] - seg_th[i])
        return float(seg_t[i] + frac * (seg_t[i + 1] - seg_t[i]))

    rising = step > 0
    t10 = first_crossing(level_10, rising)
    t90 = first_crossing(level_90, rising)

    transition_time_us = float("nan")
    if t10 is None or t90 is None:
        print("\nERROR: could not find both 10% and 90% crossings -- "
              "check the plot, the transition may not be well isolated "
              "in this capture.")
    else:
        transition_time_us = (t90 - t10) * 1e6
        print(f"\n10% crossing at t = {t10*1e3:.4f} ms")
        print(f"90% crossing at t = {t90*1e3:.4f} ms")
        print(f"Transition time (10-90%) : {transition_time_us:.2f} us")

    os_n = int(args.overshoot_window_us * 1e-6 / dt)
    seg_os = theta_p[i_trans: min(len(theta_p), i_trans + os_n)]
    if step > 0:
        overshoot_rad = float(np.max(seg_os) - after_level)
    else:
        overshoot_rad = float(after_level - np.min(seg_os))
    overshoot_deg = float(np.degrees(max(overshoot_rad, 0.0)))
    print(f"Overshoot                 : {overshoot_deg:.2f} deg")

    ripple_local = float((env_p.max() - env_p.min()) / env_p.mean() * 100)
    print(f"Envelope ripple (whole window): {ripple_local:.2f} %")

    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(t_p * 1e3, np.degrees(theta_p), color="tab:blue", lw=1.2,
            label="Measured phase")
    ax.axhline(np.degrees(before_level), color="gray", ls=":", lw=1)
    ax.axhline(np.degrees(after_level), color="gray", ls=":", lw=1)
    if t10 is not None:
        ax.axvline(t10 * 1e3, color="tab:green", ls="--", lw=1, label="10%")
    if t90 is not None:
        ax.axvline(t90 * 1e3, color="tab:red", ls="--", lw=1, label="90%")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Phase (deg)")
    ax.set_title(f"Transition detail: {transition_time_us:.1f} us "
                 f"(10-90%), overshoot {overshoot_deg:.1f} deg")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig("m3_transition_detail.png", dpi=200)
    print("\nFigure saved: m3_transition_detail.png -- check it before "
          "trusting the numbers above.")


if __name__ == "__main__":
    main()
