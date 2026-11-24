#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_nrz_rc_figures.py
======================

Builds the three QPSK NRZ-vs-raised-cosine comparison figures of
Chapter 8 (Figures 8.5, 8.6 and 8.7):

  nrz_rc_comparison_bars.png     EVM / MER / phase-error bar chart
  nrz_rc_spectrum_comparison.png NRZ and RC spectra side by side
  nrz_rc_eye_comparison.png      NRZ and RC eye diagrams side by side

It does not analyse any CSV file itself. It reads the outputs that
qpsk_analysis.py has already written:

  qpsk_nrz_results.txt, qpsk_rc_results.txt   (numbers for the bar chart)
  qpsk_nrz_spectrum.png, qpsk_rc_spectrum.png (spectrum panels)
  qpsk_nrz_eye.png, qpsk_rc_eye.png           (eye-diagram panels)

By default it looks in SDR/reproduced_plots/ (where the procedure tells
you to write your own runs) and falls back to the stored copies in
SDR/results_from_python/ and SDR/result_plots/. Output goes to
SDR/reproduced_plots/.

Usage (from the SDR folder):
  python python/make_nrz_rc_figures.py
"""

import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SDR = os.path.dirname(HERE)
OUT = os.path.join(SDR, "reproduced_plots")

INSTRUMENT_FLOOR_PCT = 0.6   # 8-bit TBS1102C EVM floor used in Chapter 8


def find(name, *folders):
    """Return the first existing path for `name` among `folders`."""
    for folder in folders:
        path = os.path.join(SDR, folder, name)
        if os.path.isfile(path):
            return path
    sys.exit(f"Cannot find {name} in {', '.join(folders)}. "
             "Run qpsk_analysis.py first (Section 5 of the procedure).")


def read_results(path):
    """Pull the four numbers needed for the bar chart from a results file."""
    text = open(path, encoding="utf-8").read()
    def grab(label):
        m = re.search(label + r"\s*:\s*([-+0-9.]+)", text)
        if not m:
            sys.exit(f"'{label}' not found in {path}")
        return float(m.group(1))
    return {
        "evm_rms": grab(r"EVM \(RMS\)"),
        "evm_peak": grab(r"EVM \(peak\)"),
        "mer": grab(r"MER"),
        "phase_rms": grab(r"Phase error, RMS"),
    }


def bar_chart(nrz, rc):
    labels = ["NRZ", "Raised-cosine"]
    dark, light = "#1f4e79", "#7fa8d9"
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 4.0))

    ax = axes[0]
    x = np.arange(2)
    w = 0.32
    b1 = ax.bar(x - w / 2, [nrz["evm_rms"], rc["evm_rms"]], w,
                label="EVM (RMS)", color=dark)
    b2 = ax.bar(x + w / 2, [nrz["evm_peak"], rc["evm_peak"]], w,
                label="EVM (peak)", color=light)
    ax.axhline(INSTRUMENT_FLOOR_PCT, color="#c00000", ls="--", lw=1.0,
               label=f"Instrument floor (~{INSTRUMENT_FLOOR_PCT}%)")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("EVM (%)")
    ax.set_title("EVM")
    ax.legend(fontsize=7.5)
    for bars in (b1, b2):
        for r in bars:
            ax.annotate(f"{r.get_height():.2f}",
                        (r.get_x() + r.get_width() / 2, r.get_height()),
                        ha="center", va="bottom", fontsize=8)

    for ax, key, title, ylabel, fmt in (
            (axes[1], "mer", "MER", "MER (dB)", "{:.2f}"),
            (axes[2], "phase_rms", "Phase Error", "Phase error, RMS (deg)", "{:.3f}")):
        bars = ax.bar(labels, [nrz[key], rc[key]], color=[dark, light], width=0.5)
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        for r in bars:
            ax.annotate(fmt.format(r.get_height()),
                        (r.get_x() + r.get_width() / 2, r.get_height()),
                        ha="center", va="bottom", fontsize=8)

    for ax in axes:
        ax.grid(True, axis="y", alpha=0.3)

    fig.suptitle("QPSK Mode 0: NRZ vs Raised-Cosine Pulse Shaping - "
                 "Software Receiver Results", fontsize=11.5)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = os.path.join(OUT, "nrz_rc_comparison_bars.png")
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print("  " + out)


def side_by_side(path_a, path_b, title, out_name):
    img_a, img_b = mpimg.imread(path_a), mpimg.imread(path_b)
    fig, axes = plt.subplots(1, 2, figsize=(14, 4.6))
    for ax, img in ((axes[0], img_a), (axes[1], img_b)):
        ax.imshow(img)
        ax.axis("off")
    fig.suptitle(title, fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0, 0.07, 1, 0.95])
    for ax, label in ((axes[0], "(a) NRZ"), (axes[1], "(b) Raised-cosine")):
        box = ax.get_position()
        fig.text(box.x0 + box.width / 2, 0.03, label, ha="center",
                 fontsize=11, color="#1f4e79", fontweight="bold")
    out = os.path.join(OUT, out_name)
    fig.savefig(out, dpi=200, facecolor="white")
    plt.close(fig)
    print("  " + out)


def main():
    os.makedirs(OUT, exist_ok=True)
    res = ("reproduced_plots", "results_from_python")
    png = ("reproduced_plots", "result_plots")

    nrz_txt = find("qpsk_nrz_results.txt", *res)
    rc_txt = find("qpsk_rc_results.txt", *res)
    print("Numbers read from:")
    print("  " + nrz_txt)
    print("  " + rc_txt)
    nrz, rc = read_results(nrz_txt), read_results(rc_txt)
    print(f"  NRZ: EVM {nrz['evm_rms']:.3f} %, MER {nrz['mer']:.2f} dB, "
          f"phase {nrz['phase_rms']:.3f} deg")
    print(f"  RC : EVM {rc['evm_rms']:.3f} %, MER {rc['mer']:.2f} dB, "
          f"phase {rc['phase_rms']:.3f} deg")

    print("Writing figures")
    bar_chart(nrz, rc)
    side_by_side(find("qpsk_nrz_spectrum.png", *png),
                 find("qpsk_rc_spectrum.png", *png),
                 "QPSK Mode 0 Spectrum: NRZ vs Raised-Cosine",
                 "nrz_rc_spectrum_comparison.png")
    side_by_side(find("qpsk_nrz_eye.png", *png),
                 find("qpsk_rc_eye.png", *png),
                 "QPSK Mode 0 Eye Diagram: NRZ vs Raised-Cosine",
                 "nrz_rc_eye_comparison.png")


if __name__ == "__main__":
    main()
