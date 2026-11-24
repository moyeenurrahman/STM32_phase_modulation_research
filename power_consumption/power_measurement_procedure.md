# Power Consumption Measurement: Procedure and Reproduction Guide

Part of the thesis *"Implementation and Performance Analysis of a Baseband Phase Modulation (PM) System with STM32 for Low-Power Wireless Communication"* (BSc EEE).

This guide explains how the **power consumption** of the STM32 baseband PM generator was measured. It tells you which firmware was flashed for each mode, how the board was wired for each of the three measurement planes (A, B and C), where the multimeter probes were placed, which Python script made each plot, and what each plot shows.

You can follow it to **repeat the measurement on your own Nucleo-F303RE**, and you can **redraw every plot in this folder without any hardware** (Section 4).

---

## 1. What is measured, and the main results

Power was measured at **three places** (called *planes*), because each plane answers a different question:

| Plane | What it measures | How | Reading |
|---|---|---|---|
| **A** | The whole Nucleo board | USB power meter in the USB line | 5.08 to 5.09 V, 0.05 A, **254 mW** (the same for every mode) |
| **B** | The STM32 chip only | 4.6 Ω shunt in the IDD jumper (JP6) | **43.6 mW** in the implemented system (mode P1) |
| **C** | The low-pass filter (op-amp) only | 458 Ω shunt in the op-amp supply line | **296 µA, 1.17 mW** |

Plane A cannot see the difference between modes. Its meter reads in steps of 10 mA, which is about 51 mW at 5.09 V, and that is larger than every difference between modes. This is why Plane B was needed.

### Headline results (Plane B, 72 MHz, mean of three runs)

| Mode (old name → new name) | What the firmware does | Current | Power |
|---|---|---|---|
| M0 → **P0** | Baseline. No DAC, DMA or timer. CPU sleeps. | 11.74 mA | **37.0 mW** |
| M1 → **P1** | PM from a precomputed buffer. CPU sleeps. **This is the implemented system.** | 14.19 mA | **43.6 mW** |
| M2 → **P2** | Same as P1, but the CPU runs an empty busy loop | 25.37 mA | **76.5 mW** |
| M3 → **P3** | PM recomputed at run time with `sinf`/`cosf` (built at `-O0`) | 30.27 mA | **87.0 mW** |

| Derived quantity | Value |
|---|---|
| Incremental power of the PM function, ΔP = P1 − P0 | **6.6 mW** |
| Saving from using sleep mode, P2 − P1 | **32.9 mW** (43.0 % of P2) |
| Saving from precomputing the waveform, P3 − P1 | **43.4 mW** (49.9 % of P3) |
| Energy per bit at 500 bit/s, Eb = ΔP / Rb | **13.2 µJ/bit** |
| Whole PM system, MCU + filter, 43.6 + 1.17 | **44.8 mW** |
| MCU share of the board total (43.6 of 254 mW) | **17.2 %** |
| Clock-scaling slope (four PLL points) | **153.7 µA/MHz**, intercept 3.145 mA, R² = 0.99986 |

### How to read these results honestly

1. **P3 is an upper bound.** It was built at `-O0` (no optimisation), and the timing margin of the callbacks was not measured. An optimised build would use less power.
2. **The multimeter is a low-cost DT830D**, not laboratory calibrated. The shunt resistance limits Plane B to about **±5.4 %** in absolute value. This error is the same in every mode, so it **cancels in the differences** (ΔP, the two savings, and Eb). The *differences* are therefore more trustworthy than the absolute values.
3. **Plane A is the board, not the chip.** The ST-LINK debugger, LEDs and regulator losses are included. The MCU is only 17.2 % of it. Do not quote 254 mW as the power of the PM system.
4. **The filter is powered from a separate battery.** Its 1.17 mW is added to the MCU power to get the system figure. It is not part of the 254 mW or the 43.6 mW.

---

## 2. Folder map

```
power_consumption/
├── power_measurement_procedure.md   <- this guide
├── firmware/
│   ├── P0_mode.c        <- baseline (was M0)
│   ├── P1_mode.c        <- implemented PM system (was M1)
│   ├── P2_mode.txt      <- busy-loop mode (was M2). Saved as .txt, copy it into main.c
│   └── P3_mode.c        <- run-time computation mode (was M3)
├── python/
│   ├── clock_scaling_plot.py   <- makes clock_scaling.png
│   └── power_figures.py        <- makes fig2, fig3 and fig4
└── result_plots/
    ├── clock_scaling.png
    ├── fig2_power_by_mode.png
    ├── fig3_energy_per_bit.png
    └── fig4_plane_breakdown.png
```

**There are no CSV files and no oscilloscope images in this folder.** Power was read from a USB meter and a multimeter, so the raw data are a few numbers per run. **All of those raw readings are written out in Section 11**, and the Python scripts contain the averaged values. The oscilloscope was used only to check that the correct waveform was present before each reading.

**The four clock-sweep firmware files (M4a to M4d) are not stored here.** They are the P1 firmware with four small changes, which are listed in full in Section 6.3.

---

## 3. Mode names: M0 to M3 are the same as P0 to P3

The four modes were first called **M0, M1, M2, M3** and were later renamed **P0, P1, P2, P3**. They are the **same modes**.

| Old name | New name | Firmware file |
|---|---|---|
| M0 | P0 | `P0_mode.c` |
| M1 | P1 | `P1_mode.c` |
| M2 | P2 | `P2_mode.txt` |
| M3 | P3 | `P3_mode.c` |

The clock-sweep modes keep their old names: **M4a, M4b, M4c, M4d** (48, 32, 16 and 8 MHz). They are the P1 firmware at lower clock speeds.

Some comments in the firmware files and some labels in the Python scripts still use the old M names (details in Section 6.3 and Section 13). This is only a naming difference. The numbers are not affected.

---

## 4. Quick reproduction of the plots (no hardware)

**Requirements:** Python 3.9 or newer with `numpy` and `matplotlib`.

```bash
pip install numpy matplotlib
cd power_consumption/python
python clock_scaling_plot.py
python power_figures.py
```

The scripts need **no input files**. The measured averages are typed into the top of each script. The images are saved in the folder you run the script from, as `.png` (300 dpi) and `.pdf`.

You should see this printed:

```
Slope  k  = 153.7 uA/MHz
Intercept = 3.145 mA
R^2       = 0.99986
Predicted at 8 MHz with PLL = 4.37 mA
Measured  at 8 MHz, PLL off = 3.59 mA
PLL current cost            = 0.78 mA
...
Incremental power dP        = 6.6 mW
Sleep-mode saving           = 32.9 mW  (43.0 %)
Precomputation saving       = 43.4 mW  (49.9 %)
MCU as fraction of board    = 17.2 %
```

---

## 5. Hardware

### 5.1 Signal chain

```
STM32 Nucleo-F303RE
  TIM6 TRGO --> DAC1_OUT1 (pin PA4), circular DMA, 2000-sample buffer
                    |
              Sallen-Key 2nd-order Butterworth low-pass filter
              R1 = R2 = 5.6 kOhm, C1 = 2 nF, C2 = 1 nF (fc about 20.1 kHz), MCP6002 op-amp
                    |
              Oscilloscope CH1 (only to check the waveform)
```

The filter has its **own battery supply** and takes no current from the Nucleo 3.3 V rail. Its only links to the board are the signal from PA4 and a common ground.

### 5.2 Fixed signal parameters

| Parameter | Value |
|---|---|
| DAC sample rate Fs | 100 kSa/s (TIM6 trigger) |
| Carrier fc | 5 kHz (20 samples per cycle) |
| Message fm | 500 Hz (200 samples per cycle) |
| Phase sensitivity kp | 1.0 rad |
| Waveform | s(t) = cos(2π·5000·t + 1.0·sin(2π·500·t)) |
| DAC mid-level, amplitude | 2048 and 1000 codes (about 1.6 V peak-to-peak) |
| Buffer length N | 2000 samples = 100 carrier cycles = 10 message cycles |

N must be a multiple of lcm(20, 200) = 200. If it is not, the waveform jumps in phase each time the DMA buffer wraps around.

### 5.3 Low-pass filter: pin-to-pin connections (MCP6002, DIP-8)

Op-amp pin numbers, looking at the chip from the top with the notch at the top-left:

| Pin | Name | Connect to |
|---|---|---|
| 1 | OUT A | Filter output. Goes to the scope CH1. Also tied to pin 2 and to one end of C1. |
| 2 | IN A− | Pin 1 (unity-gain follower) |
| 3 | IN A+ | Node B (the junction of R2 and C2) |
| 4 | VSS | Ground |
| 5 | IN B+ | Ground (a fixed DC level, so the unused second op-amp inputs do not float) |
| 6 | IN B− | Pin 7 |
| 7 | OUT B | Pin 6 (second op-amp left as a follower, output not used) |
| 8 | VDD | Battery positive, through the 458 Ω shunt (see Plane C) |

Filter network:

```
PA4 ──[ R1 5.6k ]── node A ──[ R2 5.6k ]── node B ──> pin 3 (IN A+)
                      |                       |
                      C1 = 2 nF               C2 = 1 nF
                      |                       |
                      └──> pin 1 (OUT)        └──> GND

pin 2 ──> pin 1        (unity gain)
pin 5 ──> GND          pin 6 ──> pin 7
pin 4 ──> GND          100 nF between pin 8 and pin 4 (decoupling)
```

Check: fc = 1 / (2π √(R1·R2·C1·C2)) = 20.1 kHz, and Q = √(C1/C2) / 2 = 0.707, which is the Butterworth value. **C1 (2 nF) goes from node A to the output. C2 (1 nF) goes from node B to ground.** Swapping them ruins the response.

The 100 nF decoupling capacitor sits **on the op-amp side of the shunt resistor** (between pin 8 and pin 4), so it does not disturb the voltage reading across the shunt.

---

## 6. Firmware

### 6.1 How to build one mode

1. Use an STM32CubeIDE project for the **Nucleo-F303RE**.
2. In the `.ioc` file set up:
   - **DAC1, OUT1** on pin PA4, trigger = **Timer 6 trigger out event**.
   - **DMA** for DAC1_CH1: memory to peripheral, **circular** mode, half-word data width, memory address increment on.
   - **TIM6**: Prescaler 0, Counter period (ARR) 719, trigger output = **Update Event**.
   - **Clock**: HSE in bypass mode (Section 6.2).
3. Copy the chosen firmware file into `Core/Src/main.c`.
4. Build, flash, then **follow the run rules in Section 7**.

### 6.2 Clock setup (all modes P0 to P3)

| Item | Setting |
|---|---|
| Source | **HSE 8 MHz, bypass mode** (the clock comes from the ST-LINK on the Nucleo board) |
| PLL | PREDIV /1, multiplier ×9 = 72 MHz |
| HCLK, PCLK2 | 72 MHz |
| PCLK1 (APB1 ÷2) | 36 MHz, timer clock 72 MHz |
| Flash latency | 2 wait states |
| TIM6 | PSC = 0, ARR = 719, so Fs = 72 MHz / 720 = 100.000 kSa/s |

Why HSE and not HSI: the internal HSI cannot reach 72 MHz on this chip. The PLL input from HSI is 4 MHz (HSI/2), and the PLL multiplier stops at ×16, which gives 64 MHz at most.

In the F3 HAL, the PLL divider is set with `RCC_OscInitStruct.PLL.PREDIV = RCC_PREDIV_DIV1;`.

### 6.3 The four firmware modes (Planes A and B)

| Mode | File | DAC, DMA, TIM6 | Waveform buffer | Main loop | Scope shows |
|---|---|---|---|---|---|
| **P0** (M0) | `P0_mode.c` | **Not initialised** (only `MX_GPIO_Init`) | none | `__WFI()` (sleep) | Flat line |
| **P1** (M1) | `P1_mode.c` | Running | Computed **once** at start (`Generate_Analog_PM_Buffer`) | `__WFI()` (sleep) | PM wave, 5 kHz |
| **P2** (M2) | `P2_mode.txt` | Running | Same as P1 | **Empty `while(1){}`** (CPU awake) | PM wave, same as P1 |
| **P3** (M3) | `P3_mode.c` | Running | Recomputed **all the time** in the DMA half-transfer and transfer-complete callbacks (`sinf` and `cosf` per sample) | Empty `while(1){}` | PM wave, same as P1 |

What each comparison tells you:

- **P1 − P0** = the cost of running the DAC, DMA and TIM6 (ΔP = 6.6 mW).
- **P2 − P1** = what sleep mode saves. Only the main loop differs. The waveform is identical.
- **P3 − P1** = what precomputing the waveform saves. The output is identical, because the recomputed samples have the same values as the precomputed ones.

Extra notes:

- **P3 needs the DMA interrupt enabled** in the NVIC (DMA1 channel of DAC1_CH1). Without it the callbacks never run.
- In **P3**, the buffer is filled once in `main()` before the DAC starts. After that, the half-transfer callback refills samples 0 to 999 and the transfer-complete callback refills samples 1000 to 1999.
- **Do not call `__HAL_DBGMCU_ENABLE_DBGSLEEP`.** It keeps debug clocks alive in sleep and adds current.
- The comments in `P0_mode.c`, `P2_mode.txt` and `P3_mode.c` still say "identical to M1" or "identical in M1, M2 and M3". M1 is the old name of P1 (Section 3). The **code** is correct: P2 and P3 have an empty `while(1){}` and no `__WFI()`.
- The header of `P3_mode.c` says "See the note in the source below on build optimisation". That note is not in the file. The build used for the measurement was **`-O0` (no optimisation)**, which is why P3 is treated as an upper bound.

### 6.4 Clock sweep modes (Plane B only): M4a to M4d

The firmware is **P1 with only the clock changed**. Nothing else changes: same signal, same buffer, same `__WFI()` sleep, same init order. TIM6's ARR is changed so that Fs stays at exactly 100 kSa/s. **APB1 is kept at ÷2 for every point**, so no prescaler change enters the sweep.

| Mode | HCLK | PLL multiplier | PCLK1 (APB1 ÷2) | TIM6 clock | TIM6 ARR | Flash latency |
|---|---|---|---|---|---|---|
| P1 (72 MHz) | 72 MHz | ×9 | 36 MHz | 72 MHz | 719 | 2 |
| M4a | 48 MHz | ×6 | 24 MHz | 48 MHz | 479 | 1 |
| M4b | 32 MHz | ×4 | 16 MHz | 32 MHz | 319 | 1 |
| M4c | 16 MHz | ×2 | 8 MHz | 16 MHz | 159 | 0 |
| M4d | 8 MHz | **PLL off**, HSE is the system clock | 4 MHz | 8 MHz | 79 | 0 |

To make a variant, edit only these in `SystemClock_Config()` of the P1 file: `PLL.PLLMUL`, the flash latency argument of `HAL_RCC_ClockConfig`, and (for M4d) set `PLL.PLLState = RCC_PLL_NONE` and `SYSCLKSource = RCC_SYSCLKSOURCE_HSE`. Then set the TIM6 ARR in CubeIDE (`htim6.Init.Period`). Check the "APB1 timer clocks" box in the Clock Configuration view after each change. Too few flash wait states causes a hard fault. Too many wastes power and spoils the sweep.

---

## 7. Instrument and rules for every reading

**Instruments**

- **Plane A:** USB power meter in the USB line, reads volts and amps (current resolution 10 mA).
- **Planes B and C:** DT830D digital multimeter, 3½ digits (1999 counts). DC-voltage accuracy: ±(0.5 % + 2 digits). All Plane B and C readings are **DC voltage** readings (the current is calculated from the shunt voltage), so the meter's current ranges are never used.
- Oscilloscope on the filter output (pin 1) to confirm the waveform.

**Rules for every run**

1. Flash the code from STM32CubeIDE. **Terminate the debug session** (the red square).
2. **Unplug and replug the USB cable** (power-cycle). The firmware then runs **standalone from flash**, not under the debugger. Readings taken while the debugger is attached are not valid.
3. Check the waveform on the scope. **P0 must show a flat line**. P1, P2 and P3 must show the PM wave.
4. **Warm up 5 minutes** before the first reading of a mode.
5. Take **three independent runs per mode**: power down, unplug, replug, wait 30 s, read again.
6. Write down all readings of the run before moving on.

---

## 8. Plane A: board level (USB meter)

**Connection**

```
PC / USB port ──> [ USB power meter ] ──> Nucleo USB connector (ST-LINK)
```

The IDD jumper **JP6 stays fitted**. Everything on the board is included in the reading.

**Procedure for each mode (P0, P1, P2, P3)**

1. Flash the mode, terminate the debug session, power-cycle, and check the waveform (Section 7).
2. Read **V** and **I** on the USB meter. Do this three times (three runs).
3. Power = V × I.

**Result:** 5.08 to 5.09 V and 0.05 A in **every run of every mode**, so P = **254 mW** every time. The 10 mA resolution hides all mode differences. Plane A is used for the **board-level total only**.

---

## 9. Plane B: the chip only (IDD jumper and 4.6 Ω shunt)

This is the plane used for all mode comparisons and the clock sweep.

**Idea:** the Nucleo has a jumper **JP6 (IDD)** in series with the MCU supply. Remove the jumper and put a small resistor in its place. The voltage across the resistor gives the MCU current.

**The shunt:** a 4.7 Ω resistor. Its measured resistance with the meter leads was 4.9 Ω, and the lead resistance was 0.3 Ω. So **R = 4.9 − 0.3 = 4.6 Ω** (uncertainty about ±5.4 %).

**Connections (on a breadboard, with male-to-female jumper wires)**

```
Nucleo JP6 pin 1        breadboard                           Nucleo JP6 pin 2
(supply side)           row P ──[ 4.6 Ω shunt ]── row Q      (MCU side)
      │                   │                         │              │
      └──── wire ─────────┘                         └──── wire ────┘
```

1. **Remove the JP6 jumper** from the board.
2. JP6 pin 1 (the supply side) → breadboard row **P**.
3. Put the 4.6 Ω shunt between row **P** and row **Q**.
4. Row **Q** → JP6 pin 2 (the MCU side).

**Where the multimeter probes go (three readings per run, one after another)**

| Reading | Red probe | Black probe | Range | What it gives |
|---|---|---|---|---|
| **V_shunt** | row P | row Q | 200 mV | Current: I = V_shunt / 4.6 Ω |
| **V_DD** | row Q | Nucleo GND | 20 V | The voltage the MCU really gets. Used for P = V_DD × I |
| **V_supply** | row P | Nucleo GND | 20 V | Supply before the shunt (cross-check: V_supply − V_DD should equal V_shunt) |

**Procedure for each mode**

1. Remove JP6, build the connection above.
2. Flash the mode (P0, P1, P2, P3 or M4a to M4d). Terminate the debug session. Power-cycle.
3. Check the scope waveform. Warm up 5 minutes.
4. Read V_shunt, V_DD and V_supply.
5. Repeat for three independent runs (Section 7).
6. Calculate I and P (Section 12).

**Note on scatter:** in mode P1, runs 2 and 3 gave the same shunt voltage (63.9 mV) but V_DD differed by 0.11 V (3.02 V and 3.13 V). V_supply − V_DD was 0.07 V and 0.06 V in these two runs, so the shunt reading is consistent. The V_DD difference is most likely due to variable contact resistance in the JP6 jumper-wire connections, upstream of the shunt. If a run looks odd, re-seat the connectors and repeat it.

---

## 10. Plane C: the low-pass filter only (battery and 458 Ω shunt)

**Idea:** the filter has its own supply. Put a shunt in that supply line and measure the filter current.

**Supply:** a **3.7 V Li-ion cell**, not a switching (buck) converter. A switching converter could inject ripple into the filter supply, and the same filter is used for the SNR and THD measurements. The battery also keeps the filter current separate from the Nucleo. The cell read **4.11 V** open-circuit at the start and **4.08 V** at the end.

**Shunt:** **458 Ω** (measured with the same meter before the test). The current is about 0.3 mA, so this gives a good reading of about 136 mV.

**Connections**

```
Battery (+) ── row X ──[ 458 Ω shunt ]── row Y ──> MCP6002 pin 8 (VDD)
Battery (−) ───────────────────────────────────> common GND (also Nucleo GND)
                                          100 nF between pin 8 and pin 4
```

The filter input is PA4, the filter output is pin 1, and the pin-to-pin wiring is in Section 5.3.

**Where the multimeter probes go (three readings per run)**

| Reading | Red probe | Black probe | Range | What it gives |
|---|---|---|---|---|
| **V_shunt** | row X | row Y | 200 mV | Current: I = V_shunt / 458 Ω |
| **V_pin8** | row Y | GND | 20 V | The voltage the op-amp really gets |
| **V_batt,loaded** | row X | GND | 20 V | Battery voltage under load |

Also measure the **battery open-circuit voltage** at the start and at the end of the test, to show it did not sag.

**Procedure**

1. Measure the shunt resistance first (458 Ω).
2. Build the connection above. Common ground between battery, filter and Nucleo.
3. Run firmware **P1** (the filter sees the PM waveform from PA4). Terminate the debug session, power-cycle, and check the filtered PM wave on the scope at pin 1.
4. Take the three readings for each of three runs.
5. Calculate I and P (Section 12).

The filter current is dominated by the quiescent current of the two MCP6002 amplifiers (the datasheet typical value is about 100 µA per amplifier).

---

## 11. Raw readings

### 11.1 Plane B, 72 MHz modes (R = 4.6 Ω)

| Mode | Run | V_shunt (mV) | V_DD (V) | V_supply (V) |
|---|---|---|---|---|
| M0 / P0 | 1 | 53.8 | 3.15 | 3.21 |
| | 2 | 54.5 | 3.14 | 3.19 |
| | 3 | 53.7 | 3.17 | 3.23 |
| M1 / P1 | 1 | 68.0 | 3.07 | 3.14 |
| | 2 | 63.9 | 3.02 | 3.09 |
| | 3 | 63.9 | 3.13 | 3.19 |
| M2 / P2 | 1 | 114.7 | 2.95 | 3.10 |
| | 2 | 118.5 | 3.07 | 3.19 |
| | 3 | 117.0 | 3.02 | 3.14 |
| M3 / P3 | 1 | 140.0 | 2.87 | 2.99 |
| | 2 | 139.0 | 2.92 | 3.04 |
| | 3 | 138.7 | 2.83 | 2.97 |

### 11.2 Plane B, clock sweep (R = 4.6 Ω, M1 firmware at lower clocks)

| Mode | Run | V_shunt (mV) | V_DD (V) | V_supply (V) |
|---|---|---|---|---|
| M4a (48 MHz) | 1 | 46.8 | 3.14 | 3.19 |
| | 2 | 48.5 | 3.14 | 3.19 |
| | 3 | 50.7 | 3.14 | 3.19 |
| M4b (32 MHz) | 1 | 36.6 | 3.18 | 3.22 |
| | 2 | 36.9 | 3.19 | 3.23 |
| | 3 | 37.2 | 3.18 | 3.22 |
| M4c (16 MHz) | 1 | 25.4 | 3.19 | 3.22 |
| | 2 | 25.8 | 3.20 | 3.23 |
| | 3 | 26.1 | 3.20 | 3.23 |
| M4d (8 MHz, PLL off) | 1 | 16.2 | 3.20 | 3.22 |
| | 2 | 16.5 | 3.20 | 3.22 |
| | 3 | 16.8 | 3.20 | 3.22 |

The 72 MHz point of the sweep is mode P1 (Section 11.1).

### 11.3 Plane C (R = 458 Ω, mode P1)

| Run | V_shunt (mV) | V_pin8 (V) | V_batt,loaded (V) |
|---|---|---|---|
| 1 | 135.8 | 3.95 | 4.08 |
| 2 | 135.7 | 3.95 | 4.08 |
| 3 | 135.6 | 3.94 | 4.07 |

Battery open-circuit: 4.11 V at the start, 4.08 V at the end.

### 11.4 Plane A

All modes (P0 to P3), all runs: 5.08 to 5.09 V and 0.05 A. P0 and P1 read 5.09 V, P2 read 5.08 or 5.09 V, and P3 read 5.08 V.

---

## 12. Calculations

**Plane B, per run**

```
I = V_shunt / 4.6 Ω
P = V_DD × I
```

Then average the three runs. The error bar in the plots is **half of the spread** (the largest run minus the smallest run, divided by 2).

| Mode | Mean I (mA) | Mean P (mW) | Spread of P (mW) |
|---|---|---|---|
| P0 | 11.74 | 37.0 | 0.4 |
| P1 | 14.19 | 43.6 | 3.4 (the script uses the rounded value 3.5) |
| P2 | 25.37 | 76.5 | 5.5 |
| P3 | 30.27 | 87.0 | 2.9 |

| Sweep point | Mean I (mA) | Spread of I (mA) |
|---|---|---|
| 72 MHz (P1) | 14.19 | 0.89 |
| 48 MHz (M4a) | 10.58 | 0.85 |
| 32 MHz (M4b) | 8.02 | 0.13 |
| 16 MHz (M4c) | 5.61 | 0.15 |
| 8 MHz (M4d) | 3.59 | 0.13 |

**Plane C**

```
I = V_shunt / 458 Ω = 135.7 mV / 458 Ω = 296.3 µA
P = V_pin8 × I      = 3.95 V × 296.3 µA ≈ 1.17 mW
```

**Energy per bit**

```
ΔP = P1 − P0 = 43.6 − 37.0 = 6.6 mW
Eb = ΔP / Rb
```

| Rb (bit/s) | Eb (µJ/bit) |
|---|---|
| 200 | 33.0 |
| **500 (design rate)** | **13.2** |
| 1000 | 6.6 |
| 2000 | 3.3 |

Only the **incremental** power ΔP is charged to the modulation, because the baseline 37.0 mW is the chip being powered and clocked whether it transmits or not.

**Clock-scaling fit.** A straight line is fitted by least squares through the **four PLL points only** (72, 48, 32 and 16 MHz):

```
I = I0 + k × f_HCLK        k = 153.7 µA/MHz,  I0 = 3.145 mA,  R² = 0.99986
```

k is the current that grows with clock speed (the dynamic part). I0 is the current that does not depend on the clock (the static part and the peripherals' own current). The 8 MHz point runs with the PLL off, so it is kept **off the fit**. The line predicts 4.37 mA at 8 MHz, and 3.59 mA was measured. The difference, **0.78 mA**, is attributed to switching the PLL off.

**Uncertainty.** The random scatter is the run-to-run spread above. The shunt value (±5.4 % in Plane B) and the meter accuracy are **systematic** and are the same for all modes, so they cancel in differences such as ΔP, the sleep saving and the precomputation saving.

---

## 13. Python scripts and plots

Both scripts read **no files**. The averages in Section 12 are typed into them.

| Script | Makes | Output files |
|---|---|---|
| `python/clock_scaling_plot.py` | Clock-scaling plot | `clock_scaling.png`, `.pdf` |
| `python/power_figures.py` | Three figures (Fig. 2, 3, 4) | `fig2_power_by_mode`, `fig3_energy_per_bit`, `fig4_plane_breakdown` (`.png`, `.pdf`) |

### What each plot shows

| Plot | What it shows | What to read from it |
|---|---|---|
| **`clock_scaling.png`** | MCU current (mA) against HCLK frequency (MHz). The four PLL points have error bars and a fitted line. The 8 MHz point (PLL off) is a separate red square, with the 0.78 mA PLL cost marked. | Current rises in a straight line with clock speed: **153.7 µA per MHz**. The 8 MHz point lies below the line because the PLL is off. |
| **`fig2_power_by_mode.png`** | Bar chart of MCU power for P0 to P3 (37.0, 43.6, 76.5, 87.0 mW) with error bars. Brackets mark the two savings. | Sleep saves **32.9 mW (43.0 %)**. Precomputing saves **43.4 mW (49.9 %)**. P1 is the implemented system. |
| **`fig3_energy_per_bit.png`** | Energy per bit against symbol rate (log scale) at 200, 500, 1000 and 2000 bit/s. | Eb = ΔP / Rb. At the design rate of 500 bit/s, Eb = **13.2 µJ/bit**. A faster symbol rate lowers the energy per bit. |
| **`fig4_plane_breakdown.png`** | Left: power at Plane A, B and C (log scale, 254, 43.6 and 1.17 mW). Right: the 254 mW board total split into MCU and platform overhead. | The MCU is only **17.2 %** of the board total. The ST-LINK, LEDs and regulator losses make up the other 82.8 %. The filter is small compared with the MCU. |

**Note on old labels in the scripts:** the stored `fig2_power_by_mode.png` carries the new labels P0 to P3. In `power_figures.py` the `modes` list near the top (the line starting `modes  = [`) still reads `M0` to `M3`, and the title of Figure 4 still says "M1 condition". The title of `clock_scaling_plot.py` still says "M1 firmware". Change these labels to P0 to P3 and P1 if you want the scripts to give exactly the same pictures as the stored images. The numbers are the same.

---

## 14. Limits of this measurement

1. **The multimeter is not calibrated**, and the shunt gives about ±5.4 % in absolute value. Quote **differences**, not absolute values, when comparing modes.
2. **P3 is an upper bound** (`-O0` build, callback timing margin not measured).
3. **Plane A has a 10 mA step.** It cannot separate the modes and is used only for the board total.
4. **The 8 MHz point differs in two ways from the other four**: no PLL, and the system clock comes straight from HSE. It is therefore kept off the fit.
5. **Plane C uses a battery, not the Nucleo supply**, so the filter current does not add to the Plane A or B readings. The MCU-plus-filter total of 44.8 mW is the sum of two separately measured planes.
6. **Each reading is DC (averaged by the meter).** No time-resolved current trace was recorded, so peaks (for example, the wake-ups in sleep mode) are not visible.
