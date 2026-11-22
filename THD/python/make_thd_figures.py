#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
 EXTRA THD FIGURES  --  figures 2, 3 and 5 for the thesis
===============================================================================

 Run thd_analysis.py first, then run this.  It reuses the loader and the
 spectrum code from thd_analysis.py, so the two scripts cannot drift apart.

   Figure 2  instrument floor    : function generator overlaid on T2
   Figure 3  op-amp evidence     : T8A (Point A) overlaid on T8B (Point B)
   Figure 5  THD versus amplitude, both test points, with the measured
             instrument floor drawn as a horizontal line

 Output goes to ./figures as PDF and PNG at 300 dpi.

===============================================================================
"""

import os
import glob
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from thd_analysis import (CASES, DATA_DIR, FIG_DIR, N_HARMONICS,
                          load_scope_csv, spectrum, lobe_power, unique_files,
                          analyse_run, aggregate, db)

# Which cases appear in the amplitude sweep, and which case is the floor
SWEEP_B = ["T5B", "T6B", "T7", "T2", "T8B"]
SWEEP_A = ["T5A", "T6A", "T1", "T8A"]
FLOOR_CASE = "FG"          # the function generator capture
OPERATING_POINT = 1200     # DAC codes used for the PM experiments


# =============================================================================
#  HELPERS
# =============================================================================

def case_files(name):
    """Return the de-duplicated file list for one case."""
    cfg = CASES.get(name)
    if cfg is None:
        return []
    files = sorted(glob.glob(os.path.join(DATA_DIR, cfg["pattern"])))
    files, _ = unique_files(files)
    return files


def mean_spectrum_dbc(name):
    """
    Median power spectrum of a case, in dB relative to its own fundamental.

    Each run is normalised by its OWN fundamental before the runs are
    combined.  That keeps the vertical axis meaningful even though the runs
    were taken at slightly different amplitudes, and it means two cases
    captured at different volts/div can still be compared on one axis.

    The median across runs is used only for display.  Every THD number in
    the thesis comes from the per-run analysis in thd_analysis.py, which
    detects the fundamental separately in each file.
    """
    cfg = CASES[name]
    files = case_files(name)
    if not files:
        return None, None

    stack, faxis = [], None
    for p in files:
        t, v = load_scope_csv(p)
        dt = float(np.mean(np.diff(t)))
        f, P = spectrum(v, dt)
        df = f[1] - f[0]

        band = np.abs(f - cfg["f0"]) <= max(50.0, 3 * df)
        k0 = np.where(band)[0][int(np.argmax(P[band]))]
        P1 = lobe_power(P, df, f[k0])
        if P1 <= 0:
            continue

        if faxis is None:
            faxis = f
        if len(P) != len(faxis):
            continue
        stack.append(P / P1)

    if not stack:
        return None, None

    med = np.median(np.vstack(stack), axis=0)
    return faxis, 10 * np.log10(med + 1e-30)


def case_f0(name):
    """Detected fundamental of a case, averaged over its runs."""
    cfg = CASES[name]
    files = case_files(name)
    if not files:
        return cfg["f0"]
    f0s = []
    for p in files[:4]:
        try:
            f0s.append(analyse_run(p, cfg["f0"], None, False)["f0"])
        except Exception:
            pass
    return float(np.mean(f0s)) if f0s else cfg["f0"]


def case_thd(name):
    """THD in dB and its 1-sigma spread in dB, for one case."""
    cfg = CASES.get(name)
    files = case_files(name)
    if not files:
        return None, None
    runs = []
    for p in files:
        try:
            runs.append(analyse_run(p, cfg["f0"], None, cfg["lpf_corr"]))
        except Exception:
            pass
    agg = aggregate(runs)
    if agg is None:
        return None, None
    return db(agg["thd_f"]), agg["thd_db_std"]


def save(fig, stem):
    os.makedirs(FIG_DIR, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG_DIR, f"{stem}.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {FIG_DIR}/{stem}.pdf and .png")


def mark_harmonics(ax, n=N_HARMONICS):
    """
    Mark the harmonic positions on an axis whose x units are multiples of
    the fundamental.  Labels sit at the bottom so they cannot collide with
    the legend.
    """
    for h in range(1, n + 1):
        ax.axvline(h, ls="--", alpha=0.30, lw=0.7, color="0.4")
        ax.text(h, 0.02, f"H{h}", ha="center", va="bottom", fontsize=7.5,
                color="0.3", transform=ax.get_xaxis_transform())


# =============================================================================
#  FIGURE 2  --  instrument floor
# =============================================================================

def figure2():
    print("\nFigure 2 - function generator versus STM32")
    f_fg, P_fg = mean_spectrum_dbc(FLOOR_CASE)
    f_t2, P_t2 = mean_spectrum_dbc("T2")

    if P_fg is None:
        print("  SKIPPED: no function generator files found.")
        return
    if P_t2 is None:
        print("  SKIPPED: no T2 files found.")
        return

    n_fg = len(case_files(FLOOR_CASE))
    n_t2 = len(case_files("T2"))
    f0 = case_f0("T2")

    f0_fg = case_f0(FLOOR_CASE)

    # The x axis is frequency divided by each trace's OWN fundamental.
    # Both traces then line up on H1, H2, H3 and so on even if the
    # generator was not set to exactly the same frequency as the STM32.
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(f_t2 / f0, P_t2, lw=0.8, color="#1f77b4",
            label=f"STM32 DAC, filter output (T2, {n_t2} runs, "
                  f"$f_0$ = {f0:.0f} Hz)")
    ax.plot(f_fg / f0_fg, P_fg, lw=0.8, color="#d62728", alpha=0.85,
            label=f"function generator direct to scope ({n_fg} run"
                  f"{'s' if n_fg != 1 else ''}, $f_0$ = {f0_fg:.0f} Hz)")

    mark_harmonics(ax)
    ax.set_xlabel("Frequency, as a multiple of the fundamental")
    ax.set_ylabel("Magnitude (dB relative to fundamental)")
    ax.set_title("Measured spectrum of a known low-distortion source and of "
                 "the STM32 output\nthrough the same oscilloscope chain, "
                 "$f_0$ = 1 kHz")
    ax.set_xlim(0, N_HARMONICS + 1)
    ax.set_ylim(-100, 8)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8.5, loc="upper right", framealpha=0.95)
    save(fig, "fig2_instrument_floor")


# =============================================================================
#  FIGURE 3  --  op-amp evidence at 1500 codes
# =============================================================================

def figure3():
    print("\nFigure 3 - Point A versus Point B at 1500 codes")
    f_a, P_a = mean_spectrum_dbc("T8A")
    f_b, P_b = mean_spectrum_dbc("T8B")

    if P_a is None or P_b is None:
        print("  SKIPPED: T8A or T8B files not found.")
        return

    n_a = len(case_files("T8A"))
    n_b = len(case_files("T8B"))
    f0 = case_f0("T8B")

    f0_a = case_f0("T8A")

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(f_a / f0_a, P_a, lw=0.8, color="#2ca02c",
            label=f"Point A, raw DAC output ({n_a} runs)")
    ax.plot(f_b / f0, P_b, lw=0.8, color="#d62728", alpha=0.9,
            label=f"Point B, filter output ({n_b} runs)")

    mark_harmonics(ax)
    ax.set_xlabel("Frequency, as a multiple of the fundamental")
    ax.set_ylabel("Magnitude (dB relative to fundamental)")
    ax.set_title("Effect of the MCP6002 output stage at 1500 codes, "
                 "$f_0$ = 1 kHz\nidentical amplitude at both test points")
    ax.set_xlim(0, N_HARMONICS + 1)
    ax.set_ylim(-100, 8)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8.5, loc="upper right", framealpha=0.95)
    save(fig, "fig3_opamp_pointA_vs_pointB")


# =============================================================================
#  FIGURE 5  --  THD versus amplitude, with the instrument floor
# =============================================================================

def figure5():
    print("\nFigure 5 - THD versus amplitude with instrument floor")

    series = []
    for grp, lab, mk, col in (
            (SWEEP_B, "Point B, filter output", "o", "#d62728"),
            (SWEEP_A, "Point A, raw DAC output", "s", "#2ca02c")):
        xs, ys, es = [], [], []
        for name in grp:
            if name not in CASES or not case_files(name):
                continue
            thd, std = case_thd(name)
            if thd is None:
                continue
            xs.append(CASES[name]["codes"])
            ys.append(thd)
            es.append(std)
        if xs:
            order = np.argsort(xs)
            series.append((np.array(xs)[order], np.array(ys)[order],
                           np.array(es)[order], lab, mk, col))

    if not series:
        print("  SKIPPED: no sweep data found.")
        return

    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    for xs, ys, es, lab, mk, col in series:
        ax.errorbar(xs, ys, yerr=es, marker=mk, capsize=3, lw=1.4,
                    ms=6, color=col, label=lab)

    # ---- set the vertical range from the data, with room for labels -------
    all_y = np.concatenate([y for _, y, _, _, _, _ in series])
    floor_db, floor_std = case_thd(FLOOR_CASE)
    lo = float(np.min(all_y)) - 4.0
    hi = float(np.max(all_y)) + 4.0
    if floor_db is not None:
        lo = min(lo, floor_db - 4.0)
        hi = max(hi, floor_db + 4.0)
    ax.set_ylim(lo, hi)

    # ---- the instrument floor ---------------------------------------------
    if floor_db is not None:
        ax.axhline(floor_db, ls="--", lw=1.5, color="0.25")
        # shade only the strip between the floor and the top of the axes
        ax.axhspan(floor_db, hi, color="0.5", alpha=0.10, zorder=0)
        ax.text(0.985, floor_db + 0.35,
                f"instrument floor  {floor_db:.1f} dB "
                f"(shaded: not resolvable)",
                transform=ax.get_yaxis_transform(),
                ha="right", va="bottom", fontsize=8.5, color="0.25")
    else:
        print("  note: no function generator capture, floor line omitted.")

    ax.axvline(OPERATING_POINT, ls=":", lw=1.2, color="0.4")
    ax.text(OPERATING_POINT, hi - 0.6, " PM operating point",
            ha="left", va="top", fontsize=8.5, color="0.4")

    ax.set_xlabel("DAC amplitude (codes)")
    ax.set_ylabel("THD (dB relative to fundamental)\nlower is better")
    ax.set_title("Total harmonic distortion versus DAC amplitude, "
                 "$f_0$ = 1 kHz\nharmonics 2 to 7, error bars are one "
                 "standard deviation over runs")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9, loc="lower left", framealpha=0.95)

    save(fig, "fig5_thd_vs_amplitude_with_floor")


# =============================================================================

def main():
    print("=" * 70)
    print(" EXTRA THD FIGURES")
    print("=" * 70)
    figure2()
    figure3()
    figure5()
    print("\nDone.")


if __name__ == "__main__":
    main()
