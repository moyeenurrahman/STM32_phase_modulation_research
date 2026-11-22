"""
make_summary_figures.py  -  Chapter 7 summary figures for the SNR chapter.

Unlike snr_method_a.py and snr_method_b_sinad.py, this script does NOT read
CSV files. Figures 4, 6 and 10 are CROSS-CASE comparisons - they summarise
numbers that already came out of six separate analysis runs, so the numbers
are entered once below rather than recomputed.

Produces three files:
    fig4_inband_snr_summary.png       - in-band SNR across N1-N6
    fig6_sinad_enob_improvement.png   - SINAD and ENOB, Point A vs Point B
    fig10_power_config_comparison.png - power configuration A vs B

============================================================================
BEFORE USING: check every number in the RESULTS block below against your
own analysis output. The values here are pre-filled from the results
already produced earlier in the SNR measurement session, but you should
verify each one against your saved script output before trusting the
figures - a transcription error here would silently produce a wrong plot.
============================================================================
"""

import numpy as np
import matplotlib.pyplot as plt

# =============================================================
# RESULTS  -  edit these to match your own analysis output
# =============================================================

# ---- Figure 4: in-band SNR across all six cases ----
# "standard" = plain in-band SNR from snr_method_a.py
# "robust"   = the robust variant (outlier bins excluded)
# band       = signal band half-width used, for the modulated cases, shown
#              as an annotation since the value is not comparable without it
CASES = [
    # label,           standard,  robust
    ("N2\nPoint A\n(tone)",   29.11, 39.82),
    ("N3\nPoint B\n(tone)",   38.62, 41.09),
    ("N4\nPM kp=1.0",         32.18, 38.82),
    ("N5\nPM kp=2.4048",      30.45, 38.30),
    ("N6\nBPSK",              18.71, 35.86),
]

# Oscilloscope's own SNR limit in this configuration (Sample mode, 8-bit).
# Printed by snr_method_a.py as "Scope's own SNR limit". Values were close
# to 44.8 dB across all cases; the mean is used as a single reference line.
SCOPE_LIMIT_DB = 44.8

# Predicted DAC in-band SNR (ideal-converter calculation, NOT attainable
# with an 8-bit oscilloscope). Shown as an annotation only, since plotting
# it to scale would need an axis break at these levels.
PREDICTED_A_DB = 69.3
PREDICTED_B_DB = 72.9


# ---- Figure 6: SINAD / ENOB improvement across the filter ----
SINAD_A, SINAD_B = 20.69, 32.32
ENOB_A,  ENOB_B  = 3.92,  5.85


# ---- Figure 10: power configuration comparison ----
# Config C (laptop on mains) was assessed qualitatively only - oscillation
# was directly visible in the time domain, so no capture set was taken.
# It is shown as a labelled gap rather than a bar.
POWER_CONFIGS = [
    # label,                          in_band_snr, robust, outlier_bins
    ("A\ncharger,\nlaptop battery",   38.62, 41.09, 11),
    ("B\nlaptop USB,\nlaptop battery",38.57, 41.21, 9),
]
CONFIG_C_NOTE = "C: laptop USB,\nlaptop mains\n\u2014 oscillation\nobserved,\nnot quantified"

SAVE_DPI = 300


# =============================================================
# FIGURE 4 - in-band SNR summary across N1-N6
# =============================================================

def make_figure4():
    labels   = [c[0] for c in CASES]
    standard = [c[1] for c in CASES]
    robust   = [c[2] for c in CASES]

    x = np.arange(len(labels))
    w = 0.36

    fig, ax = plt.subplots(figsize=(10, 6))

    b1 = ax.bar(x - w/2, standard, w, label="In-band SNR (standard)",
               color="#4C72B0")
    b2 = ax.bar(x + w/2, robust, w, label="In-band SNR (robust, "
               "outliers excluded)", color="#DD8452")

    for bars in (b1, b2):
        for rect in bars:
            h = rect.get_height()
            ax.annotate(f"{h:.1f}", xy=(rect.get_x() + rect.get_width()/2, h),
                        xytext=(0, 3), textcoords="offset points",
                        ha="center", va="bottom", fontsize=8)

    ax.axhline(SCOPE_LIMIT_DB, color="tab:red", linestyle="--",
              linewidth=1.2, alpha=0.7)
    ax.text(len(labels) - 0.55, SCOPE_LIMIT_DB + 0.6,
           f"Oscilloscope resolution limit \u2248 {SCOPE_LIMIT_DB:.1f} dB "
           "(Sample mode, 8-bit)",
           color="tab:red", fontsize=8, ha="right")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("In-band SNR (dB)")
    ax.set_title("In-band SNR across all measured conditions\n"
                "N4\u2013N6 use a modulation-dependent band (not IEEE 1241 "
                "SNR) \u2014 see Table 3")
    ax.set_ylim(0, max(SCOPE_LIMIT_DB, max(robust)) + 8)
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, axis="y", alpha=0.3)

    fig.text(0.01, 0.01,
            f"Predicted DAC in-band SNR (ideal converter, not attainable "
            f"with this instrument): Point A \u2248 {PREDICTED_A_DB:.1f} dB, "
            f"Point B \u2248 {PREDICTED_B_DB:.1f} dB.",
            fontsize=7, color="dimgray", ha="left")

    plt.tight_layout(rect=[0, 0.03, 1, 1])
    plt.savefig("fig4_inband_snr_summary.png", dpi=SAVE_DPI)
    print("Saved fig4_inband_snr_summary.png")
    plt.close(fig)


# =============================================================
# FIGURE 6 - SINAD / ENOB improvement, Point A vs Point B
# =============================================================

def make_figure6():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 5))

    # --- SINAD ---
    bars = ax1.bar(["Point A\n(raw DAC)", "Point B\n(filtered)"],
                   [SINAD_A, SINAD_B], color=["#4C72B0", "#55A868"],
                   width=0.55)
    for rect, v in zip(bars, [SINAD_A, SINAD_B]):
        ax1.annotate(f"{v:.2f} dB", xy=(rect.get_x() + rect.get_width()/2, v),
                    xytext=(0, 4), textcoords="offset points",
                    ha="center", fontsize=9)
    delta_sinad = SINAD_B - SINAD_A
    arrow_y1 = max(SINAD_A, SINAD_B) + 3.5
    ax1.annotate("", xy=(1, arrow_y1), xytext=(0, arrow_y1),
                arrowprops=dict(arrowstyle="->", color="dimgray",
                                connectionstyle="arc3,rad=-0.3"))
    ax1.text(0.5, arrow_y1 + 2.0, f"+{delta_sinad:.1f} dB",
            ha="center", fontsize=9, color="dimgray")
    ax1.set_ylabel("SINAD (dB)")
    ax1.set_title("SINAD")
    ax1.set_ylim(0, arrow_y1 + 6)
    ax1.grid(True, axis="y", alpha=0.3)

    # --- ENOB ---
    bars2 = ax2.bar(["Point A\n(raw DAC)", "Point B\n(filtered)"],
                    [ENOB_A, ENOB_B], color=["#4C72B0", "#55A868"],
                    width=0.55)
    for rect, v in zip(bars2, [ENOB_A, ENOB_B]):
        ax2.annotate(f"{v:.2f} bits", xy=(rect.get_x() + rect.get_width()/2, v),
                    xytext=(0, 4), textcoords="offset points",
                    ha="center", fontsize=9)
    delta_enob = ENOB_B - ENOB_A
    arrow_y2 = max(ENOB_A, ENOB_B) + 0.7
    ax2.annotate("", xy=(1, arrow_y2), xytext=(0, arrow_y2),
                arrowprops=dict(arrowstyle="->", color="dimgray",
                                connectionstyle="arc3,rad=-0.3"))
    ax2.text(0.5, arrow_y2 + 0.4, f"+{delta_enob:.1f} bits",
            ha="center", fontsize=9, color="dimgray")
    ax2.set_ylabel("ENOB (bits)")
    ax2.set_title("Effective Number of Bits")
    ax2.set_ylim(0, arrow_y2 + 1.2)
    ax2.grid(True, axis="y", alpha=0.3)

    fig.suptitle("SINAD and ENOB improvement across the reconstruction "
                 "filter (case N2 \u2192 N3)", fontsize=11)
    fig.text(0.5, 0.01,
            "Improvement of this size is inconsistent with DAC quantisation "
            "noise, which passes the filter largely unattenuated below "
            "20 kHz; it indicates the residual at Point A was dominated by "
            "wideband interference above the cutoff.",
            ha="center", fontsize=7.5, color="dimgray", wrap=True)

    plt.tight_layout(rect=[0, 0.06, 1, 0.94])
    plt.savefig("fig6_sinad_enob_improvement.png", dpi=SAVE_DPI)
    print("Saved fig6_sinad_enob_improvement.png")
    plt.close(fig)


# =============================================================
# FIGURE 10 - power configuration comparison
# =============================================================

def make_figure10():
    labels        = [c[0] for c in POWER_CONFIGS]
    inband        = [c[1] for c in POWER_CONFIGS]
    robust        = [c[2] for c in POWER_CONFIGS]
    outliers      = [c[3] for c in POWER_CONFIGS]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

    x = np.arange(len(labels))
    w = 0.36

    # --- left panel: SNR ---
    b1 = ax1.bar(x - w/2, inband, w, label="Standard", color="#4C72B0")
    b2 = ax1.bar(x + w/2, robust, w, label="Robust", color="#DD8452")
    for bars in (b1, b2):
        for rect in bars:
            h = rect.get_height()
            ax1.annotate(f"{h:.2f}", xy=(rect.get_x()+rect.get_width()/2, h),
                        xytext=(0, 3), textcoords="offset points",
                        ha="center", fontsize=8)
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, fontsize=8)
    ax1.set_ylabel("In-band SNR (dB)")
    ax1.set_title("SNR by power configuration")
    ax1.set_ylim(0, max(robust) + 8)
    ax1.legend(fontsize=8, loc="upper left")
    ax1.grid(True, axis="y", alpha=0.3)

    # --- right panel: outlier bin count ---
    b3 = ax2.bar(x, outliers, 0.5, color="#C44E52")
    for rect, v in zip(b3, outliers):
        ax2.annotate(str(v), xy=(rect.get_x()+rect.get_width()/2, v),
                    xytext=(0, 3), textcoords="offset points",
                    ha="center", fontsize=9)
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels, fontsize=8)
    ax2.set_ylabel("Outlier bins discarded (count)")
    ax2.set_title("Discrete interference by power configuration")
    ax2.set_ylim(0, max(outliers) + 4)
    ax2.grid(True, axis="y", alpha=0.3)

    fig.suptitle("Case N3: comparison of STM32/laptop power configuration "
                 "(all captures on laptop battery except Config C)",
                 fontsize=10)
    fig.text(0.5, 0.01,
            "Config C (laptop USB, laptop on mains charger) produced "
            "directly observable time-domain oscillation; a capture set "
            "was judged unnecessary and is omitted from both panels.",
            ha="center", fontsize=7.5, color="dimgray", wrap=True)

    plt.tight_layout(rect=[0, 0.06, 1, 0.93])
    plt.savefig("fig10_power_config_comparison.png", dpi=SAVE_DPI)
    print("Saved fig10_power_config_comparison.png")
    plt.close(fig)


# =============================================================
# RUN
# =============================================================

if __name__ == "__main__":
    make_figure4()
    make_figure6()
    make_figure10()
    print("\nAll three figures generated. Check each against the RESULTS "
         "block at the top of this script before using in the thesis.")
