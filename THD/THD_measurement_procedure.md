# THD Measurement: Procedure and Reproduction Guide

Part of the thesis *"Implementation and Performance Analysis of a Baseband Phase Modulation (PM) System with STM32 for Low-Power Wireless Communication"* (BSc EEE, RUET).

This guide explains how the **total harmonic distortion (THD)** of the STM32 baseband signal generator was measured. It tells you which firmware was flashed, how the oscilloscope was set, which CSV files were recorded, which Python script turns each set of CSV files into a result, and what each plot and oscilloscope image shows.

You can follow it to **reproduce every number and plot in this folder from the CSV files alone, without any hardware.** Section 3 does this in about two minutes.

---

## 1. What is measured, and the main results

The STM32 makes a pure sine wave with its 12-bit DAC. The signal is measured at two places:

- **Test Point A:** the raw DAC output (a staircase), pin PA4.
- **Test Point B:** the output of a Sallen-Key low-pass filter (the clean signal).

THD is measured at two test frequencies:

| Frequency | Why |
|---|---|
| **1 kHz** (main test) | All harmonics 2 to 7 (2 to 7 kHz) lie inside the filter passband, so the filter does not change them. |
| **5 kHz** (operating carrier) | This is the carrier used in the PM experiments. Harmonics 2 to 7 (10 to 35 kHz) lie near and above the filter corner (about 20 kHz), so the filter response must be removed (de-embedded). |

**Definition used.** THD_F (IEEE Std 1241) is the root of the summed power of harmonics 2 to 7, divided by the power of the fundamental. THD_R (against the total RMS) is also printed. The two agree to four decimals here, because the harmonics are tiny.

### Headline results (reproduced from this folder, file `csv/results_thd.csv`)

| Case | Point | f0 | DAC amplitude | Runs | THD (%) | THD (dB) | Spread (1 σ, dB) | Largest harmonic |
|---|---|---|---|---|---|---|---|---|
| **T2** (headline) | B | 1 kHz | 1200 codes | 8 | **0.4759** | **−46.45** | 0.78 | H3 |
| T1 | A | 1 kHz | 1200 | 8 | 0.5613 | −45.02 | 1.14 | H3 |
| T5B | B | 1 kHz | 500 | 3 | 0.4940 | −46.13 | 0.84 | H3 |
| T6B | B | 1 kHz | 800 | 3 | 0.3343 | −49.52 | 0.28 | H2 |
| T7 | B | 1 kHz | 1000 | 3 | 0.3901 | −48.18 | 0.78 | H2 |
| T8B | B | 1 kHz | 1500 | 3 | **1.3284** | **−37.53** | 0.08 | H4 |
| T5A | A | 1 kHz | 500 | 3 | 0.5357 | −45.42 | 1.17 | H3 |
| T6A | A | 1 kHz | 800 | 3 | 0.5296 | −45.52 | 0.76 | H6 |
| T8A | A | 1 kHz | 1500 | 3 | 0.6097 | −44.30 | 0.57 | H2 |
| T3 | A | 5 kHz | 1200 | 32 | 0.8185 | −41.74 | 0.53 | H6 |
| T4 (filter de-embedded) | B | 5 kHz | 1200 | 30 | 0.8244 | −41.68 | 1.09 | H7 |
| **FG** (instrument check) | direct | 1 kHz | n/a | 1 | 0.7153 | **−42.91** | n/a | H5 |

### How to read these results honestly

1. **Most numbers are upper bounds, not exact values.** The oscilloscope is 8-bit and the function-generator check (FG) shows that the *oscilloscope chain alone* reads about **−42.9 dB** of "THD" on a source that is not the STM32. Any STM32 result near or below −43 dB is therefore limited by the instrument. The STM32 is **at least as good as** the value shown, and the true figure is lower. The script marks such cases with *"Report THD as an UPPER BOUND"*.
2. **Why the instrument is limited.** The High-Resolution mode was **not** active, so the record has only 2 to 3 distinct code levels on a quiet signal (20 mV steps at 500 mV/div, 8 mV steps at 200 mV/div). There is no dither, so the quantiser adds its own spurs. The sampling ratio (31.25 kSa/s for a 1 kHz tone) also repeats every 125 samples, which puts spurs on multiples of 250 Hz, which include every harmonic.
3. **T8B is a real result, and it is the one clearly above the floor.** At 1500 codes the positive peaks are flat at Point B only. Point A (T8A) is clean at the same amplitude. This points to the MCP6002 op-amp input stage (Section 10). The script prints a warning for it (measured amplitude 1 % below theory), which is expected because the peaks are flat.
4. **The 200 mV/div cases (T5, T6, T7) have no matching floor capture.** They use only the local noise floor beside each harmonic. Treat their small differences from T2 as within the measurement limit.
5. **At 5 kHz the filter correction matters.** Without correcting for the filter, T4 would read −45.82 dB. With the Butterworth response removed it reads **−41.68 dB** (+4.14 dB). The correction uses the *design* cutoff (20.1 kHz). The measured filter sweep was not completed, so this is an assumption (Section 10).

---

## 2. Folder map

```
THD/
├── THD_measurement_procedure.md   <- this guide
├── csv/                           <- oscilloscope records, 2000 points each
│   ├── T1/    8 files (1.csv ... 8.csv)    1 kHz, 1200 codes, Point A
│   ├── T2/    8 files                      1 kHz, 1200 codes, Point B   (headline)
│   ├── T5A/   3 files                      1 kHz,  500 codes, Point A   (200 mV/div)
│   ├── T5B/   3 files                      1 kHz,  500 codes, Point B   (200 mV/div)
│   ├── T6A/   3 files                      1 kHz,  800 codes, Point A   (200 mV/div)
│   ├── T6B/   3 files                      1 kHz,  800 codes, Point B   (200 mV/div)
│   ├── T7/    3 files                      1 kHz, 1000 codes, Point B   (200 mV/div)
│   ├── T8A/   3 files                      1 kHz, 1500 codes, Point A
│   ├── T8B/   3 files                      1 kHz, 1500 codes, Point B   (flat tops)
│   ├── T9a/   8 files                      noise floor, 4 ms/div, 500 mV/div
│   ├── T9b/   3 files                      noise floor, 4 ms/div, 500 mV/div
│   ├── T9c/  32 files                      noise floor, 400 us/div, 500 mV/div (same files as SNR/csv/N1)
│   ├── fg.csv                              function generator straight into the scope (1 file)
│   ├── results_thd.csv                     the result table written by thd_analysis.py
│   └── python thd.txt                      text printed by an earlier run of thd_analysis.py
├── firmware/
│   ├── THD_main.c.c               <- main.c for the THD builds (copy into CubeIDE as main.c)
│   └── dma.c.c                    <- dma.c to go with it (copy into CubeIDE as dma.c)
├── python/
│   ├── thd_analysis.py            <- main analysis: all THD numbers, 5 plots
│   ├── make_thd_figures.py        <- 3 extra plots (needs thd_analysis.py in the same folder)
│   └── filter_response.py         <- filter sweep fit (needs a sweep file that is NOT in this repo)
├── result_plots/                  <- the plots made by the scripts (Section 8)
└── oscilloscope_images/           <- scope screenshots (Section 9)
```

**Two cases live in the SNR folder.** The 5 kHz cases were recorded once, in the SNR session, and reused:

| THD case | Same files as | Folder |
|---|---|---|
| **T3** (5 kHz, Point A) | SNR case N2 | `SNR/csv/N2/` (32 files) |
| **T4** (5 kHz, Point B) | SNR case N3 | `SNR/csv/N3/` (32 files; 2 duplicates dropped, 30 used) |
| **T9c** (5 kHz floor) | SNR case N1 | `SNR/csv/N1/` (32 files), also copied to `THD/csv/T9c/` |

The SNR firmware (`SNR/firmware/SNR.c`) made these three signals. See `SNR/SNR_measurement_procedure.md`.

---

## 3. Quick reproduction (no hardware)

**Requirements:** Python 3.9 or newer with `numpy`, `scipy` and `matplotlib`.

```bash
pip install numpy scipy matplotlib
```

(Tested with Python 3.13, numpy 2.5, scipy 1.18, matplotlib 3.11.)

1. Clone the repository and open a terminal in `THD/python/`.
2. Run the main analysis:
   ```bash
   python thd_analysis.py
   ```
3. Run the extra figures:
   ```bash
   python make_thd_figures.py
   ```

`thd_analysis.py` is already set up for this repository layout (`DATA_DIR = "../csv"`). It reads the T3 and T4 files from `../../SNR/csv/N2` and `../../SNR/csv/N3`, and the function-generator file from `../csv/fg.csv`. It takes about a minute.

**What you should see.**

- The terminal prints one block per case, then the amplitude-sweep table.
- A new file `results_thd.csv` appears in the folder you ran from. It must be **identical** to `THD/csv/results_thd.csv`. The reference values are in Section 1.
- A new folder `figures/` appears with `fig_thd_*.png` and `.pdf` (5 plots) from the first script, and `fig2_*`, `fig3_*`, `fig5_*` (3 plots) from the second.
- The first block printed for T4 says `dropped 2 duplicate file(s): 3.csv, 6.csv`. This is expected.

Outputs are written to the folder you run from. Nothing in the repository is overwritten.

> **If your copy of `thd_analysis.py` still has `floor="T9b"` on the T5A, T5B, T6A, T6B and T7 lines, change those five to `floor=None`.** With `"T9b"` the script gives different numbers for those five cases (about −42.3, −46.4, −48.2, −42.3 and −46.4 dB, because they then sit on a floor measured at a different volts/div). The published table uses `floor=None`. This was checked by running both versions.

---

## 4. The test cases

| Case | What it is | Firmware | Probe at | Volts/div | Timebase | CSV folder | Records | External floor used |
|---|---|---|---|---|---|---|---|---|
| **T1** | Raw DAC linearity | B1 (1 kHz, 1200) | A | 500 mV | 4 ms | `csv/T1` | 8 | T9a |
| **T2** | **Headline result** | B1 (1 kHz, 1200) | B | 500 mV | 4 ms | `csv/T2` | 8 | T9a |
| **T3** | 5 kHz, raw | SNR `MODE_CARRIER` | A | 500 mV | 400 µs | `SNR/csv/N2` | 32 | T9c |
| **T4** | 5 kHz, filter removed | SNR `MODE_CARRIER` | B | 500 mV | 400 µs | `SNR/csv/N3` | 30 of 32 | T9c |
| **T5A / T5B** | Amplitude 500 codes | B2 | A / B | **200 mV** | 4 ms | `csv/T5A`, `csv/T5B` | 3 + 3 | none |
| **T6A / T6B** | Amplitude 800 codes | B3 | A / B | **200 mV** | 4 ms | `csv/T6A`, `csv/T6B` | 3 + 3 | none |
| **T7** | Amplitude 1000 codes | B4 | B | **200 mV** | 4 ms | `csv/T7` | 3 | none |
| **T8A** | Amplitude 1500 codes | B5 | A | 500 mV | 4 ms | `csv/T8A` | 3 | T9a |
| **T8B** | Amplitude 1500 codes (flat tops) | B5 | B | 500 mV | 4 ms | `csv/T8B` | 3 | T9a |
| **T9a** | Noise floor for the 1 kHz cases | B7a (DAC at mid-code) | B | 500 mV | 4 ms | `csv/T9a` | 8 | n/a |
| **T9b** | Second noise floor, 1 kHz | B7a | B | 500 mV | 4 ms | `csv/T9b` | 3 | n/a (analysed, not used) |
| **T9c** | Noise floor for the 5 kHz cases | SNR `MODE_DC_DMA` | B | 500 mV | 400 µs | `csv/T9c` | 32 | n/a |
| **FG** | Instrument distortion check | none (function generator) | direct | 500 mV | 4 ms | `csv/fg.csv` | 1 | T9a |

**Why the floor cases exist.** Without them you cannot tell whether a small harmonic belongs to the signal or to the oscilloscope. Each signal case is compared with a floor captured with the DAC at a constant level.

**Why T5, T6 and T7 are at 200 mV/div.** A smaller amplitude fills more of the screen at 200 mV/div and gives finer code steps (8 mV instead of 20 mV). The 500 mV/div floor captures (T9a, T9b) do not match this setting, so these five cases have `floor=None`.

**Why 500 mV/div for T1, T2, T8.** At 200 mV/div the 1200-code and 1500-code signals clip the screen. The 5 kHz cases (T3, T4) use 500 mV/div because the SNR measurement used it.

**Why T8 was measured at both points.** At 1500 codes the signal at Point B has flat positive peaks. Point A at the same amplitude is round. Taking both proves that the flattening comes from the filter stage, not from the DAC.

---

## 5. Hardware and firmware

### 5.1 Signal chain

```
STM32 Nucleo-F303RE
  TIM6 TRGO --> DAC1_OUT1 (PA4, output buffer on), DMA1 Channel 3, circular, 1000-sample buffer
                    |
                    +--> Test Point A  (raw DAC staircase)
                    |
              Sallen-Key 2nd-order Butterworth low-pass filter
              R1 = R2 = 5.6 kOhm, C1 = 2 nF, C2 = 1 nF (fc about 20.1 kHz), MCP6002 op-amp
              (op-amp supply measured at 4.05 V)
                    |
                    +--> Test Point B  (filtered output)
```

### 5.2 Fixed parameters

| Parameter | Value |
|---|---|
| DAC sample rate fs | 100 kSa/s (TIM6 on an 8 MHz clock, ARR = 79, prescaler 0) |
| Buffer | 1000 samples, circular DMA, half-word transfers |
| Mid-level | 2048 codes |
| Clock | HSI at 8 MHz, no PLL |
| Samples per cycle | 100 (1 kHz), 20 (5 kHz) |
| Expected peak-to-peak | 2 × codes × 3.3 V / 4095 (1.934 V at 1200 codes) |

The buffer always holds a **whole number of sine cycles** (10 cycles at 1 kHz, 50 at 5 kHz). If it did not, the waveform would jump every 10 ms when the DMA wraps around, and the jump would add harmonics that have nothing to do with the DAC. The firmware stops the build with `#error` if the numbers do not divide exactly, or if the amplitude would clip.

### 5.3 The firmware: one file, several builds

`firmware/THD_main.c.c` is a complete `main.c` for an STM32CubeIDE project (Nucleo-F303RE, with DAC1, DMA, TIM6 and GPIO generated by CubeMX). **Three lines in the configuration block select the build:**

```c
#define OUTPUT_MODE         MODE_SINE      /* MODE_SINE, MODE_DC_MIDPOINT or MODE_OFF */
#define F0_HZ               1000U          /* fundamental, Hz                         */
#define DAC_AMPLITUDE       1200U          /* peak, in DAC codes                      */
```

| Build | `OUTPUT_MODE` | `F0_HZ` | `DAC_AMPLITUDE` | Used for |
|---|---|---|---|---|
| **B1** | `MODE_SINE` | 1000 | 1200 | T1, T2 |
| **B2** | `MODE_SINE` | 1000 | 500 | T5A, T5B |
| **B3** | `MODE_SINE` | 1000 | 800 | T6A, T6B |
| **B4** | `MODE_SINE` | 1000 | 1000 | T7 |
| **B5** | `MODE_SINE` | 1000 | 1500 | T8A, T8B |
| B6 | `MODE_SINE` | 2000 | 1200 | Planned (T10), **not captured** |
| **B7a** | `MODE_DC_MIDPOINT` | n/a | n/a | T9a, T9b (DAC on, constant code 2048, no DMA) |
| B7b | `MODE_OFF` | n/a | n/a | Alternative floor (DAC and DMA off), not used in the final data |
| n/a | SNR firmware | 5000 | 1200 | T3, T4 (SNR cases N2, N3) and T9c (N1): see `SNR/firmware/SNR.c` |

**Design choices that matter for a distortion measurement:**

- The main loop sleeps in `__WFI()`, with SysTick suspended. No printing, no LED, no timer interrupt reaches the supply while the signal is measured.
- The DMA half-transfer and transfer-complete interrupts are switched off **after** `HAL_DAC_Start_DMA()` (the HAL turns them on again each time it is called). Otherwise the buffer wrap would raise an interrupt at 100 Hz, and the harmonics of 100 Hz land on exactly 1, 2 and 3 kHz.
- The sine table is rounded to the nearest code. Each build recomputes it at compile time from `F0_HZ` and `DAC_AMPLITUDE`.
- This firmware produces **a pure sine only**. THD must never be taken from a PM waveform, because the wanted sidebands cannot be told apart from harmonic distortion.

### 5.4 The `dma.c` file

`firmware/dma.c.c` replaces CubeMX's generated `dma.c`. It removes two unused memory-to-memory DMA channels and the DMA2 clock (less supply current and switching noise), and keeps only the DMA1 clock and the Channel 3 interrupt line. The DAC's own DMA handle `hdma_dac1_ch1` is created by CubeMX inside `dac.c`. Check these settings there:

| Setting | Value |
|---|---|
| Instance | `DMA1_Channel3` |
| Direction | `DMA_MEMORY_TO_PERIPH` |
| Peripheral and memory alignment | Half-word (both) |
| Peripheral increment / memory increment | Disable / Enable |
| Mode | `DMA_CIRCULAR` |
| Priority | High |

Wrong circular mode or wrong alignment silently ruins the output. Remove the memory-to-memory entries in CubeMX too (System Core, DMA), otherwise the next code generation adds them back.

### 5.5 Flashing

1. Create a CubeIDE project for the Nucleo-F303RE with DAC1 (OUT1 on PA4), TIM6 (TRGO on update event), DMA1 Channel 3 and GPIO. Keep CubeMX peripheral init as **separate files** (`dac.c`, `dma.c`, `gpio.c`, `tim.c`).
2. Copy `firmware/THD_main.c.c` over `Core/Src/main.c`, and `firmware/dma.c.c` over `Core/Src/dma.c`. The two files carry a double extension (`.c.c`) in this repository, so rename them when you copy.
3. Set `OUTPUT_MODE`, `F0_HZ` and `DAC_AMPLITUDE` for the build you need.
4. **Project, Clean**, then **Project, Build**, so the table is regenerated.
5. Flash with **Run** (not Debug), then unplug the debugger USB and power the board as for the capture (Section 6.2).

---

## 6. Capture procedure (oscilloscope)

### 6.1 Equipment

Tektronix **TBS1102C** oscilloscope (8-bit, 2000-point records), 1× passive probe, a function generator (for the FG check), a computer with TekVISA and OpenChoice Desktop (or a USB drive, if allowed), the STM32 board and the filter board.

### 6.2 Wiring and power

1. Keep the wires between the Nucleo and the filter board short.
2. PA4 goes to the filter input (this node is **Test Point A**). The filter output is **Test Point B**.
3. Use the **shortest probe ground lead**, on the shared ground.
4. Run the computer on battery during the capture, and keep chargers and other switching supplies away.

### 6.3 Oscilloscope settings

| Setting | Value | Reason |
|---|---|---|
| Channel | CH1 | n/a |
| Probe | 1× (switch **and** menu) | Full amplitude to the 8-bit converter |
| Coupling | **AC** | Removes the DC offset (about 1.65 V), so the signal fills the screen |
| Vertical scale | 500 mV/div (T1, T2, T3, T4, T8A, T8B, T9a, T9b, T9c, FG); **200 mV/div** (T5A, T5B, T6A, T6B, T7) | 500 mV/div avoids clipping at 1200 and 1500 codes |
| Acquisition | **Sample** | Not Hi-Res, not Average |
| Averaging | **Off** | Averaging would hide the noise and the run-to-run spread |
| Bandwidth limit | 20 MHz on | Cuts front-end noise |
| Trigger | Edge (rising), CH1 | n/a |
| Timebase | **4 ms/div** for 1 kHz; **400 µs/div** for 5 kHz | See the table below |

| Quantity | 1 kHz cases (4 ms/div) | 5 kHz cases (400 µs/div) |
|---|---|---|
| Sample rate | 31.25 kSa/s | 312.5 kSa/s |
| Sample interval in the CSV | 32 µs | 3.2 µs |
| Record length | 2000 points = 64 ms | 2000 points = 6.4 ms |
| Cycles of the fundamental in a record | 64 | 32 |
| FFT bin width | 15.63 Hz | 156.25 Hz |
| Nyquist | 15.6 kHz (harmonics up to H7 = 7 kHz fit) | 156 kHz (H7 = 35 kHz fits) |

Both settings put an **exact whole number of cycles** in the record (64 and 32). The scripts also read the sample interval from each CSV, so you do not need to type it in.

**Note on the CSV header.** The `Horizontal Scale` line in each CSV shows 6.4 ms (and 640 µs for T9c). This is the record length divided by 10 divisions, not the front-panel setting. The panel was at 4 ms/div (and 400 µs/div).

### 6.4 Saving records

1. Press **Single** on the scope and **wait until it shows Stop**.
2. Save the record as CSV (OpenChoice: *Waveform Data Capture, CH1, Get Data, Save as CSV*).
3. Name the files `1.csv`, `2.csv` and so on, one folder per case. **Press Single before every save.** Saving twice without a new acquisition gives two identical files. This happened in the T4 data (3.csv equals 2.csv, and 6.csv equals 4.csv), and the script removes such duplicates automatically.

### 6.5 Order of work

1. Flash **B1**. Capture **T2** (probe at B), then **T1** (probe at A). Take the T2 screenshot.
2. Flash **B2**, **B3**, **B4** in turn with the scope at 200 mV/div. Capture T5A and T5B, T6A and T6B, then T7.
3. Flash **B5** with the scope at 500 mV/div. Capture **T8A**, then **T8B**. Take both screenshots (flat tops at B only).
4. Flash **B7a**. Capture the floors **T9a** and **T9b** at Point B. Take the T9a screenshot.
5. Disconnect the STM32 and feed the function generator (1 kHz, about 1.9 Vpp, to match the STM32 amplitude) directly into the scope. Capture **FG**. Take the screenshot.
6. **T3, T4, T9c** come from the SNR session. They were not captured again.

### 6.6 CSV file format

Each file is one 2000-point record in Tektronix format. The first rows carry the settings (`Record Length`, `Sample Interval`, `Vertical Scale`, `Horizontal Scale`) in the first columns, and the data is in the last two columns:

```
Record Length,2000,Points,-32.00000152E-3,-100.0000089E-3,
Sample Interval,31.99999992E-6,s,-31.96800128E-3,80.000005662E-3,
...
```

The last two columns are time (s) and voltage (V). The script finds them by checking which column rises steadily (time) and which does not (voltage), and it reads the sample interval from the header.

---

## 7. Analysis: which script, which CSV folder, in what order

Run everything from `THD/python/` (Section 3). No setting needs editing for the repository layout.

### 7.1 Order of runs

1. `python thd_analysis.py` (needs all CSV folders). Writes `results_thd.csv` and 5 plots.
2. `python make_thd_figures.py` (needs `thd_analysis.py` in the same folder, and the same CSV folders). Writes 3 plots.
3. `filter_response.py` is **not run**: its input (a frequency sweep of the filter) was not recorded. See Section 7.4.

### 7.2 `thd_analysis.py`

The `CASES` table at the top ties each case to its CSV folder, fundamental, DAC amplitude, test point and floor:

| Case | `pattern` (relative to `THD/csv`) | `f0` | `codes` | `floor` | `lpf_corr` |
|---|---|---|---|---|---|
| T9a, T9b, T9c | `T9a/*.csv`, `T9b/*.csv`, `T9c/*.csv` | 1000 / 1000 / 5000 | n/a | n/a | no |
| FG | `fg.csv` | 1000 | n/a | T9a | no |
| T1, T2 | `T1/*.csv`, `T2/*.csv` | 1000 | 1200 | T9a | no |
| T5B, T6B, T7 | `T5B/*.csv`, `T6B/*.csv`, `T7/*.csv` | 1000 | 500, 800, 1000 | **none** | no |
| T5A, T6A | `T5A/*.csv`, `T6A/*.csv` | 1000 | 500, 800 | **none** | no |
| T8A, T8B | `T8A/*.csv`, `T8B/*.csv` | 1000 | 1500 | T9a | no |
| T3 | `../../SNR/csv/N2/*.csv` | 5000 | 1200 | T9c | no |
| T4 | `../../SNR/csv/N3/*.csv` | 5000 | 1200 | T9c | **yes** |

Other settings: harmonics 2 to 7 (`N_HARMONICS = 7`), main-lobe half width 4 bins (`LOBE_BINS = 4`), floor margin 6 dB (`FLOOR_MARGIN_DB = 6.0`), filter cutoff `LPF_FCUT = 20100.0` Hz, `VDDA = 3.3`.

**What the script does, step by step:**

1. **Floors first.** It reads T9a, T9b and T9c and measures the floor power at each harmonic frequency.
2. **Load and clean.** It reads every CSV in a case and drops byte-identical duplicates (MD5 check).
3. **Spectrum.** Each 2000-point record is windowed with a **Blackman-Harris** window (side lobes far below the harmonics, which Hann would not give) and transformed with an FFT.
4. **Find the fundamental** in each record separately, within ±50 Hz of the nominal value. The STM32 and scope clocks are independent, so it is never exactly 1000.00 Hz in every record.
5. **Harmonic power.** The power of the fundamental and of each harmonic is the **sum over the whole main lobe** (±4 bins), not one bin.
6. **Floor test.** The floor beside each harmonic is estimated locally (median of nearby bins). If an external floor is set for the case, the higher of the two is used. A harmonic that is **less than 6 dB above the floor** is not trusted. Its floor value is used instead, so it counts as an **upper bound**. The printout shows `measured` or `at/below floor` and how many runs passed.
7. **Filter de-embedding (T4 only).** Each harmonic and the fundamental are divided by the Butterworth gain at their own frequency, so the result describes the signal before the filter.
8. **Combine the runs.** The **power ratios** of the runs are averaged (not the spectra, which would smear the harmonics). THD_F is the square root of the summed ratios. The spread is the standard deviation of the per-run values in dB.
9. **Checks printed with every case:** detected fundamental, expected and measured Vpp (a warning appears if the measured amplitude is more than 1 % below theory), and the **crest factor** (peak divided by RMS). A pure sine gives 1.414. A flattened peak lowers it.
10. **Extra comparison.** T3 and T4 are also analysed on an **8-run subset**, so that they can be compared fairly with the 1 kHz cases, which have 8 runs. Result: T3 −41.38 dB, T4 −41.45 dB.

**Reference printout.** `csv/python thd.txt` holds the text from an earlier run. Its THD values match `results_thd.csv`. (Its first lines describe T9b as "200 mV/div". This is a wrong label in the script's note: T9b was recorded at 500 mV/div, as the CSV headers show.)

### 7.3 `make_thd_figures.py`

Uses the same loader, the same spectrum code and the same `CASES` as `thd_analysis.py`, so the two cannot drift apart. It produces:

| Figure | Cases used | Purpose |
|---|---|---|
| `fig2_instrument_floor` | FG and T2 | Shows the scope chain alone against the STM32 |
| `fig3_opamp_pointA_vs_pointB` | T8A and T8B | Shows the op-amp flat-top effect |
| `fig5_thd_vs_amplitude_with_floor` | T5A to T8B | THD against amplitude, with the −42.9 dB floor line from FG |

In these plots the x-axis is frequency divided by each trace's own fundamental, and the vertical axis is dB relative to the fundamental.

### 7.4 `filter_response.py` (not used for the published numbers)

It fits a cutoff frequency, Q and passband gain to a filter sweep (`filter_sweep.csv` with columns `freq_Hz, vin_pp, vout_pp`) and prints the filter correction factors. The sweep was **not recorded**, so this script has no input and no output in this repository. The T4 correction uses the **design** Butterworth response (`LPF_FCUT = 20100` Hz, 2nd order) in `thd_analysis.py`. If a sweep is added later, rerun this script and compare the fitted fc with 20.1 kHz. A large difference would change the T4 number.

---

## 8. Result plots (`result_plots/`)

The repository file names drop the `fig_` or `figN_` prefix that the scripts add.

| Repository file | Script output name | Script | What it shows | What to look for |
|---|---|---|---|---|
| `thd_spectrum_T2.png` | `fig_thd_spectrum_T2` | `thd_analysis.py` | Spectrum of the headline case (Point B, 1 kHz, 1200 codes), 0 to 8 kHz, with H1 to H7 marked and the T9a noise floor overlaid | A sharp peak at 1 kHz, small peaks at H3 and H5 (about −50 to −52 dBc) that rise above the floor, and a floor of about −70 to −80 dBc |
| `thd_A_vs_B.png` | `fig_thd_A_vs_B` | `thd_analysis.py` | T1 (Point A) and T2 (Point B) spectra overlaid, 0 to 10 kHz | The two traces nearly match. Point A also shows two lines near 5.6 and 7.6 kHz that are about −40 dBc. They are **not harmonics** (they do not sit on multiples of 1 kHz), and the filter removes them at Point B. |
| `thd_harmonics_T2.png` | `fig_thd_harmonics_T2` | `thd_analysis.py` | Bar chart of H2 to H7 for T2, in dBc | Hatched bars (H2, H4, H6) are at or below the noise floor, so only upper bounds. H3, H5 and H7 are measured. |
| `thd_vs_amplitude.png` | `fig_thd_vs_amplitude` | `thd_analysis.py` | THD (dB) against DAC amplitude at Point A and Point B, error bars = 1 σ over runs | Point B stays between −46 and −50 dB up to 1200 codes, then jumps to −37.5 dB at 1500 codes. Point A stays near −45 dB throughout. |
| `thd_vs_amplitude_with_floor.png` | `fig5_thd_vs_amplitude_with_floor` | `make_thd_figures.py` | The same plot, with the instrument floor (−42.9 dB, from FG) as a dashed line and the region **below** it shaded | Every point except T8B lies **inside the shaded region**, which is the formal reason they are upper bounds (a THD lower than the instrument's own distortion cannot be resolved). T8B, above the line, is the only point clearly outside the instrument limit. |
| `instrument_floor.png` | `fig2_instrument_floor` | `make_thd_figures.py` | Function generator (red) and STM32 Point B (blue) spectra on the same axes | The function generator alone shows H2, H3 and H5 at about −50 dBc. So the scope chain adds harmonics of its own, and it is not more linear than the STM32. |
| `opamp_pointA_vs_pointB.png` | `fig3_opamp_pointA_vs_pointB` | `make_thd_figures.py` | T8A (green) and T8B (red) at the same 1500-code amplitude | Point B has strong H3 to H7 lines (H4 at about −45 dBc). Point A is much cleaner at those harmonics. This is the evidence for the op-amp cause. |
| `thd_5kHz_correction.png` | `fig_thd_5kHz_correction` | `thd_analysis.py` | T4 harmonic levels (H2 to H7 at 10 to 35 kHz) before (raw) and after removing the filter response | The correction raises H4 to H7, because the filter attenuates them. THD changes from −45.82 dB to −41.68 dB. |

---

## 9. Oscilloscope images (`oscilloscope_images/`)

Screenshots from the TBS1102C. They are for illustration. **The analysis uses the CSV files only.** The screenshots were taken at **different timebase settings** from the CSV records (so that a few cycles fill the screen), so the sample rate shown on the screen does not match the CSV.

| File | Case | Screen settings | What it shows |
|---|---|---|---|
| `T2_cope.bmp` | T2 (Point B, 1200 codes) | 500 mV/div, 400 µs/div | A clean, round sine. Pk-pk 1.960 V, frequency about 1.005 kHz, RMS 680.5 mV. |
| `T8A_scope.png` | T8A (Point A, 1500 codes) | 500 mV/div, 200 µs/div | Both peaks round (staircase visible). Pk-pk 2.460 V, max 2.890 V, min 430 mV. |
| `T8B_scope.png` | T8B (Point B, 1500 codes) | 500 mV/div, 200 µs/div | **Flat positive peaks**, rounded negative peaks. Pk-pk 2.380 V, max 2.830 V, min 450 mV. |
| `T9a_scope.png` | T9a (floor) | 500 mV/div, 400 µs/div | An almost flat line with a few LSB steps. Pk-pk 40 mV. Shows the oscilloscope noise floor. (The "frequency" reading of 3.255 kHz is meaningless on noise.) |
| `fg.png` | FG | 500 mV/div, 1 ms/div | A function-generator 1 kHz sine fed straight into the scope. Pk-pk 1.820 V. It is the reference for the −42.9 dB instrument floor. |

---

## 10. Findings behind the numbers

**Instrument floor (FG).** A function generator into the scope gives THD_F = −42.91 dB (0.7153 %), with H2, H3 and H5 clearly present. So the measurement cannot show an STM32 THD better than about −43 dB. T1, T2 and the sweep points sit near or below that line, and are therefore **upper bounds**.

**Filter and amplitude (T2 against T1).** T2 (Point B) is 1.4 dB better than T1 (Point A). The difference is of the same size as the run-to-run spread (0.8 to 1.1 dB), so it is **not a firm result**.

**Op-amp flattening at 1500 codes (T8B).** At 1500 codes the output swing is 0.44 to 2.86 V. The MCP6002 is powered from 4.05 V. A rail-to-rail op-amp with a complementary input pair changes from one input pair to the other near Vdd − 1.2 V (2.85 V here), and a 2.86 V peak reaches that point. The positive peaks flatten (the screenshot shows the flat tops, and the crest factor at the positive side is 1.406 against 1.419 on the negative side). Point A at the same amplitude is round, and its THD stays near −44 dB. So the distortion in T8B comes from the filter stage, and it is the one effect clearly bigger than the instrument limit. **Conclusion for the PM experiments:** keep the amplitude at or below 1200 codes (T2, about 0.68 to 2.62 V), where this does not happen.
**5 kHz (T3 and T4).** T3 (−41.74 dB) and T4 after correction (−41.68 dB) agree. Both are within about 1.2 dB of the 1 kHz instrument floor (−42.9 dB), and the script flags most of their harmonics as at or below the noise floor, so treat them as upper bounds too.

**Cross-check with SINAD (SNR folder).** For the 5 kHz signal, SINAD should equal the combination of SNR and THD:

`SINAD = −10 log10( 10^(−SNR/10) + 10^(THD/10) )`

Using the wideband SNR (0.1 to 148 kHz) from `SNR/` and THD from this folder:

| Point | Wideband SNR | THD (this folder) | Predicted SINAD | SINAD from sine fit (SNR folder) |
|---|---|---|---|---|
| A (T3 with N2) | 20.76 dB | −41.74 dB | 20.73 dB | 20.69 dB |
| B (T4 with N3) | 32.40 dB | −41.68 dB | 31.92 dB | 32.32 dB |

The two independent measurements agree within 0.4 dB. At these SNR values the noise is the main limit, and THD at −42 dB adds only a small share.

---

## 11. Known limits

- **8-bit oscilloscope, no High-Resolution mode, no dither.** The −42.9 dB instrument floor (FG) is the best resolvable THD. Anything lower is an upper bound. Recording with High-Resolution mode, or with a 12-bit oscilloscope, would lower this limit.
- **Only 3 runs for T5 to T8, and 1 run for FG.** The spread values for 3 runs are rough. The 8 runs of T1, T2 and 30 to 32 runs of T3, T4 are firmer.
- **No 200 mV/div floor.** T5, T6 and T7 rely on the local noise floor only.
- **Filter response not measured.** T4 uses the design cutoff.
- **B6 (2 kHz) was planned but not captured.**
- **Quantiser spurs.** At 1 kHz the sampling ratio repeats every 125 samples, so spurs appear at multiples of 250 Hz. These can fall on harmonics. This is part of the instrument limit above.
- **The 5 kHz data comes from the SNR session** (different day settings and a different firmware file), as described in Section 2.

---

## 12. Quick checklist to reproduce everything

1. `pip install numpy scipy matplotlib`
2. `cd THD/python`
3. `python thd_analysis.py`. Check `results_thd.csv` against `../csv/results_thd.csv` (it must match exactly) and compare the plots in `figures/` with `../result_plots/`.
4. `python make_thd_figures.py`. Check `fig2_*`, `fig3_*`, `fig5_*`.
5. To redo the measurement on hardware: Section 5 (firmware builds B1 to B5, B7a) and Section 6 (oscilloscope settings and order of work).
