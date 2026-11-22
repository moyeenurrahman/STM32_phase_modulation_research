"""
Clock-scaling characteristic for the baseband PM generator
STM32 Nucleo-F303RE, Plane B (IDD jumper, 4.6 ohm shunt)

Plots MCU supply current against system clock frequency for the M1 firmware
(PM active, precomputed buffer, CPU in Sleep mode), fits a straight line to
the four PLL-based points, and shows the PLL-disabled point separately.

    I_MCU = I0 + k * f_HCLK

    k  = dynamic current coefficient (uA/MHz)
    I0 = static and peripheral current that does not scale with clock
"""

import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Measured data (Plane B, mean of three independent runs)
# ---------------------------------------------------------------------------

# Points running from the PLL
f_pll = np.array([72.0, 48.0, 32.0, 16.0])          # MHz
i_pll = np.array([14.19, 10.58, 8.02, 5.61])        # mA
e_pll = np.array([0.89, 0.85, 0.13, 0.15]) / 2.0    # half-spread, mA

# Point running directly from HSE with the PLL disabled
f_nopll = np.array([8.0])
i_nopll = np.array([3.59])
e_nopll = np.array([0.13]) / 2.0

# ---------------------------------------------------------------------------
# Least-squares fit through the PLL points only
# ---------------------------------------------------------------------------

k, i0 = np.polyfit(f_pll, i_pll, 1)

i_fit = i0 + k * f_pll
ss_res = np.sum((i_pll - i_fit) ** 2)
ss_tot = np.sum((i_pll - np.mean(i_pll)) ** 2)
r2 = 1.0 - ss_res / ss_tot

i_pred_8 = i0 + k * 8.0
pll_cost = i_pred_8 - i_nopll[0]

print(f"Slope  k  = {k*1000:.1f} uA/MHz")
print(f"Intercept = {i0:.3f} mA")
print(f"R^2       = {r2:.5f}")
print(f"Predicted at 8 MHz with PLL = {i_pred_8:.2f} mA")
print(f"Measured  at 8 MHz, PLL off = {i_nopll[0]:.2f} mA")
print(f"PLL current cost            = {pll_cost:.2f} mA")

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 11,
    "axes.linewidth": 0.8,
})

fig, ax = plt.subplots(figsize=(6.5, 4.5))

f_line = np.linspace(0, 80, 200)
ax.plot(f_line, i0 + k * f_line,
        color="0.35", linewidth=1.2, zorder=1,
        label=fr"Fit: $I = {i0:.2f} + {k*1000:.1f}\times10^{{-3}} f$   ($R^2 = {r2:.4f}$)")

ax.errorbar(f_pll, i_pll, yerr=e_pll,
            fmt="o", markersize=6, color="#1f4e79",
            markerfacecolor="#1f4e79", markeredgecolor="black",
            markeredgewidth=0.6, capsize=3, linestyle="none", zorder=3,
            label="Measured, PLL enabled")

ax.errorbar(f_nopll, i_nopll, yerr=e_nopll,
            fmt="s", markersize=7, color="#c0392b",
            markerfacecolor="#c0392b", markeredgecolor="black",
            markeredgewidth=0.6, capsize=3, linestyle="none", zorder=3,
            label="Measured, PLL disabled (HSE direct)")

ax.plot([8.0, 8.0], [i_nopll[0], i_pred_8],
        color="#c0392b", linestyle=":", linewidth=1.2, zorder=2)
ax.annotate(f"PLL cost\n{pll_cost:.2f} mA",
            xy=(8.0, (i_nopll[0] + i_pred_8) / 2),
            xytext=(17, 3.3),
            fontsize=9, color="#c0392b",
            arrowprops=dict(arrowstyle="->", color="#c0392b", linewidth=0.8))

ax.set_xlabel("System clock frequency, $f_{HCLK}$ (MHz)")
ax.set_ylabel("Microcontroller supply current, $I_{MCU}$ (mA)")
ax.set_title("Clock-scaling characteristic of the baseband PM generator\n"
             "(M1 firmware: PM active, precomputed buffer, CPU in Sleep mode)",
             fontsize=10.5, pad=12)

ax.set_xlim(0, 80)
ax.set_ylim(0, 16)
ax.grid(True, which="major", linestyle="--", linewidth=0.5, alpha=0.5)
ax.legend(loc="upper left", frameon=True, framealpha=1.0,
          edgecolor="0.7", fontsize=9)

fig.tight_layout()
fig.savefig("clock_scaling.png", dpi=300, bbox_inches="tight")
fig.savefig("clock_scaling.pdf", bbox_inches="tight")

print("\nSaved clock_scaling.png and clock_scaling.pdf")
