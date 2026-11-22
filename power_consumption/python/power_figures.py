"""
Power consumption figures for Chapter 7
STM32 Nucleo-F303RE baseband PM generator

Figure 2 : MCU power by operating mode at 72 MHz (Plane B)
Figure 3 : Energy per bit against symbol rate
Figure 4 : Power breakdown across measurement planes

Style matched to the clock-scaling figure.
"""

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 11,
    "axes.linewidth": 0.8,
})

NAVY   = "#1f4e79"
RED    = "#c0392b"
TEAL   = "#2a7f7f"
GREY   = "#7f7f7f"

# ===========================================================================
# FIGURE 2 - MCU power by operating mode, Plane B, 72 MHz
# ===========================================================================

modes  = ["M0\nBaseline", "M1\nPM, sleeping", "M2\nPM, busy loop", "M3\nRuntime cosf()"]
power  = np.array([37.0, 43.6, 76.5, 87.0])          # mW
perr   = np.array([0.4, 3.5, 5.5, 2.9]) / 2.0        # half-spread, mW
colors = [GREY, NAVY, RED, RED]

fig, ax = plt.subplots(figsize=(6.5, 4.6))

bars = ax.bar(modes, power, yerr=perr, capsize=4, width=0.62,
              color=colors, edgecolor="black", linewidth=0.7, zorder=3)

# Value labels
for i, (p, e) in enumerate(zip(power, perr)):
    ax.text(i, p + e + 2.0, f"{p:.1f} mW",
            ha="center", va="bottom", fontsize=9.5, fontweight="bold")

# Highlight the implemented system, placed inside the bar to avoid collision
ax.text(1, 43.6 / 2, "implemented\nsystem", ha="center", va="center",
        fontsize=8.5, color="white", style="italic")

# Saving annotations, drawn as brackets above the bars
def bracket(ax, x1, x2, y, text, color):
    ax.plot([x1, x1, x2, x2], [y - 1.5, y, y, y - 1.5],
            color=color, linewidth=1.0, zorder=4)
    ax.text((x1 + x2) / 2, y + 0.8, text, ha="center", va="bottom",
            fontsize=9, color=color)

bracket(ax, 1, 2, 96, "sleep-mode saving\n32.9 mW  (43.0 %)", RED)
bracket(ax, 1, 3, 111, "precomputation saving\n43.4 mW  (49.9 %)", RED)

ax.axhline(37.0, color=GREY, linestyle=":", linewidth=1.0, zorder=1)
ax.text(3.42, 37.0, "baseline", fontsize=8.5, color=GREY,
        va="center", ha="left")

ax.set_ylabel("Microcontroller power, $P_{MCU}$ (mW)")
ax.set_title("Microcontroller power by operating mode\n"
             "(Plane B, 72 MHz, mean of three independent runs)",
             fontsize=10.5, pad=12)
ax.set_ylim(0, 128)
ax.set_xlim(-0.6, 3.9)
ax.grid(True, axis="y", linestyle="--", linewidth=0.5, alpha=0.5, zorder=0)
ax.set_axisbelow(True)

fig.tight_layout()
fig.savefig("fig2_power_by_mode.png", dpi=300, bbox_inches="tight")
fig.savefig("fig2_power_by_mode.pdf", bbox_inches="tight")
plt.close(fig)
print("Figure 2 saved")

# ===========================================================================
# FIGURE 3 - Energy per bit against symbol rate
# ===========================================================================

dP = 43.6 - 37.0                       # incremental power, mW
Rb = np.array([200.0, 500.0, 1000.0, 2000.0])
Eb = dP / Rb * 1000.0                  # uJ/bit

# Continuous curve
Rb_line = np.logspace(np.log10(150), np.log10(3000), 200)
Eb_line = dP / Rb_line * 1000.0

fig, ax = plt.subplots(figsize=(6.5, 4.5))

ax.plot(Rb_line, Eb_line, color="0.45", linewidth=1.2, zorder=1,
        label=r"$E_b = \Delta P / R_b$,   $\Delta P = 6.6$ mW")

ax.plot(Rb, Eb, "o", markersize=7, color=NAVY,
        markerfacecolor=NAVY, markeredgecolor="black",
        markeredgewidth=0.6, zorder=3, linestyle="none",
        label="Reported operating points")

# Highlight the design symbol rate
ax.plot(500, dP / 500 * 1000, "o", markersize=11, markerfacecolor="none",
        markeredgecolor=RED, markeredgewidth=1.6, zorder=4)
ax.annotate("design symbol rate\n500 bit/s,  13.2 µJ/bit",
            xy=(500, 13.2), xytext=(830, 22),
            fontsize=9, color=RED, ha="left",
            arrowprops=dict(arrowstyle="->", color=RED, linewidth=0.9,
                            connectionstyle="arc3,rad=-0.2"))

for r, e in zip(Rb, Eb):
    if r != 500:
        ax.annotate(f"{e:.1f} µJ/bit", xy=(r, e), xytext=(9, 7),
                    textcoords="offset points", ha="left",
                    fontsize=8.5, color=NAVY)

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Symbol rate, $R_b$ (bit/s)")
ax.set_ylabel("Energy per bit, $E_b$ (µJ/bit)")
ax.set_title("Energy per bit against symbol rate\n"
             "(incremental power at the microcontroller plane)",
             fontsize=10.5, pad=12)
ax.set_xlim(150, 3200)
ax.set_ylim(2.2, 70)

# Explicit ticks at the reported operating points, plain number formatting
ax.set_xticks([200, 500, 1000, 2000])
ax.set_xticklabels(["200", "500", "1000", "2000"])
ax.set_xticks([], minor=True)
ax.set_yticks([3, 5, 10, 20, 40, 60])
ax.set_yticklabels(["3", "5", "10", "20", "40", "60"])
ax.set_yticks([], minor=True)

ax.grid(True, which="major", linestyle="--", linewidth=0.5, alpha=0.45)
ax.legend(loc="lower left", frameon=True, framealpha=1.0,
          edgecolor="0.7", fontsize=9)

fig.tight_layout()
fig.savefig("fig3_energy_per_bit.png", dpi=300, bbox_inches="tight")
fig.savefig("fig3_energy_per_bit.pdf", bbox_inches="tight")
plt.close(fig)
print("Figure 3 saved")

# ===========================================================================
# FIGURE 4 - Power breakdown across measurement planes
# ===========================================================================

fig, (axL, axR) = plt.subplots(1, 2, figsize=(9.0, 4.3),
                               gridspec_kw={"width_ratios": [1.0, 1.25]})

# --- Left panel: the three planes on a log scale -------------------------
plane_labels = ["Plane A\nboard level", "Plane B\nMCU only", "Plane C\nLPF stage"]
plane_power  = np.array([254.0, 43.6, 1.17])
plane_colors = [GREY, NAVY, TEAL]

axL.bar(plane_labels, plane_power, width=0.6,
        color=plane_colors, edgecolor="black", linewidth=0.7, zorder=3)

for i, p in enumerate(plane_power):
    label = f"{p:.0f} mW" if p >= 100 else (f"{p:.1f} mW" if p >= 10 else f"{p:.2f} mW")
    axL.text(i, p * 1.35, label, ha="center", va="bottom",
             fontsize=9.5, fontweight="bold")

axL.set_yscale("log")
axL.set_ylabel("Power (mW)")
axL.set_ylim(0.5, 900)
axL.set_title("Measured power at each plane\n(M1 condition, logarithmic scale)",
              fontsize=10, pad=10)
axL.grid(True, axis="y", which="both", linestyle="--",
         linewidth=0.5, alpha=0.45, zorder=0)
axL.set_axisbelow(True)

# --- Right panel: composition of the board-level total -------------------
mcu      = 43.6
lpf      = 1.17          # externally powered, shown for context only
overhead = 254.0 - 43.6

axR.barh([0.2], [mcu], color=NAVY, edgecolor="black", linewidth=0.7,
         height=0.45, zorder=3, label=f"STM32F303RE  ({mcu:.1f} mW,  17.2 %)")
axR.barh([0.2], [overhead], left=[mcu], color=GREY, edgecolor="black",
         linewidth=0.7, height=0.45, zorder=3,
         label=f"Platform overhead  ({overhead:.1f} mW,  82.8 %)")

axR.text(mcu / 2, 0.2, "MCU", ha="center", va="center",
         fontsize=9, color="white", fontweight="bold")
axR.text(mcu + overhead / 2, 0.2,
         "ST-LINK debugger, LEDs,\nregulator losses",
         ha="center", va="center", fontsize=8.5, color="white")

axR.set_xlim(0, 275)
axR.set_ylim(-1.15, 0.55)
axR.set_yticks([])
axR.set_xlabel("Power (mW)")
axR.set_title("Composition of the board-level total\n"
              "(Plane A = 254 mW)", fontsize=10, pad=10)
axR.legend(loc="lower center", bbox_to_anchor=(0.5, 0.02),
           frameon=True, framealpha=1.0, edgecolor="0.7", fontsize=8.5)
axR.grid(True, axis="x", linestyle="--", linewidth=0.5, alpha=0.45, zorder=0)
axR.set_axisbelow(True)

fig.tight_layout()
fig.savefig("fig4_plane_breakdown.png", dpi=300, bbox_inches="tight")
fig.savefig("fig4_plane_breakdown.pdf", bbox_inches="tight")
plt.close(fig)
print("Figure 4 saved")

# ===========================================================================
# Printed check values
# ===========================================================================
print()
print(f"Incremental power dP        = {dP:.1f} mW")
print(f"Sleep-mode saving           = {76.5-43.6:.1f} mW  ({(76.5-43.6)/76.5*100:.1f} %)")
print(f"Precomputation saving       = {87.0-43.6:.1f} mW  ({(87.0-43.6)/87.0*100:.1f} %)")
print(f"MCU as fraction of board    = {43.6/254.0*100:.1f} %")
print(f"LPF as fraction of MCU+LPF  = {1.17/(43.6+1.17)*100:.1f} %")
for r, e in zip(Rb, Eb):
    print(f"  Rb = {r:6.0f} bit/s  ->  Eb = {e:5.1f} uJ/bit")
