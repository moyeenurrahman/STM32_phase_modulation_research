# Output Frequency Accuracy: Measurement and Reproduction Guide

Part of the thesis *"Implementation and Performance Analysis of a Baseband Phase Modulation (PM) System with STM32 for Low-Power Wireless Communication"* (BSc EEE, RUET).

This guide explains how the **output frequency accuracy** of the STM32 baseband signal generator was measured. It tells you which firmware was flashed, how the oscilloscope was set, which CSV files were recorded, and which Python script turns each set of CSV files into a result. You can follow it to **reproduce every number and plot from the CSV files in this folder**, without any hardware.

> **Test names.** The tests were first called **T1 … T7** and later renamed **F1 … F7**. They are the same tests: **T1 = F1, T2 = F2, … T7 = F7**. The firmware files and this guide use "F". The CSV files, Python scripts, and result files still use "T". Both names appear below.

---

## 1. What is measured, and the main result

**Quantity:** the frequency error of the generated carrier, in parts per million (ppm):

`error_ppm = (f_measured − f_set) / f_set × 1 000 000`

**Key idea.** The output frequency comes from a chain of exact integer divisions (PLL, prescaler, timer reload, DMA buffer):

`f_out = f_source × (M_PLL / D_PLL) × 1 / [(PSC+1)(ARR+1)] × (n_cycles / N)`

Every factor except `f_source` is an exact integer ratio, so **the output error equals the clock-source error**. This has three consequences:

1. Testing several frequencies (1, 2, 5, 10 kHz) is a *consistency check*. All of them must show the same error in ppm, with the same sign.
2. The real result is the clock source: the internal **HSI** RC oscillator versus the **HSE** crystal reference.
3. The measurement is limited by the reference. The oscilloscope timebase is ±20 ppm and the HSE crystal is ±20 ppm. So **HSI can be measured**, but **HSE can only be bounded**.

### Headline results (all reproducible from this folder)

| Quantity | Result |
|---|---|
| HSI error, DAC output (F1, 12 captures) | **+3106.8 ppm** (+0.31 %), sample σ 87.2 ppm, pooled within-frequency σ 81.4 ppm |
| HSI error is the same at all frequencies (F1) | One-way ANOVA across 1/2/5/10 kHz: p = 0.277 (no frequency dependence) |
| HSE error, DAC output (F2, 12 captures) | mean +82.8 ppm, **all captures below 137 ppm**. This is a bound, not a value (see Section 2) |
| Ratiometric HSI − HSE at 5 kHz (F3) | **+3125.7 ppm**, which cancels the scope timebase error. Uncertainty: ±20 ppm (crystal) plus about ±60 ppm run-to-run |
| Clock pin PA8 (MCO), HSI (F4) | +3047.2 ppm, versus +3106.8 ppm through the DAC path. They agree (Welch t = −0.91, p = 0.44), so the DAC/DMA/timer add no error |
| PM mode, β = 1 rad (F6) | carrier +72.1 ppm, versus +89.2 ppm unmodulated at 5 kHz (Welch t = −0.99). Modulation does not shift the carrier |
| Filter check (F7) | filtered output +36.05 ppm; raw DAC staircase unreliable (see Section 8). The filter does not shift the frequency within ±60 ppm |
| Short-term stability over about 31 min (F5) | HSI scatter 86.1 ppm, HSE 17.9 ppm; no significant warm-up trend (HSI p = 0.275, HSE p = 0.599) |

---

## 2. The reference problem

| Element | Accuracy |
|---|---|
| STM32F303 HSI RC oscillator | ±1 % (10 000 ppm), at room temperature |
| HSE (8 MHz from the Nucleo ST-LINK, bypass mode) | ±20 ppm |
| Tektronix TBS1102C timebase | ±20 ppm |
| Zero-crossing estimator (this work) | about 1 ppm at the record length used |

- HSI error (about 3000 ppm) is 150 times larger than the scope's 20 ppm, so the scope is a valid reference for it.
- HSE error is the same size as the scope's 20 ppm. A scope alone cannot measure it, so the result is reported as a bound (every HSE capture is below 137 ppm, and every per-frequency mean is below 110 ppm).
- Comparing HSI and HSE in one setup (F3) cancels the scope timebase error, because the same timebase measures both.

---

## 3. Equipment and wiring

| Item | Use |
|---|---|
| STM32 Nucleo-F303RE | Waveform generator, DAC1_OUT1 on **PA4**, TIM6 trigger → DMA |
| Sallen-Key Butterworth low-pass filter (MCP6002 op-amp, f_c ≈ 20.1 kHz, Q = 0.707) | Reconstruction filter. **Reported measurements are taken at its output** |
| Tektronix TBS1102C oscilloscope | 2 channels, 8-bit, **fixed 2000-point record**, ±20 ppm timebase |
| 10× passive probe | The scope channel attenuation must also be set to **10×** (a mismatch gave a wrong 0.232 Vpp early on) |
| USB flash drive or OpenChoice PC software | CSV export |
| STM32CubeIDE | Build and flash |
| Python 3 | numpy ≥ 2.0 and matplotlib (see Section 9) |

**Signal path:** `PA4 (DAC) → Sallen-Key LPF → scope CH1 → CSV → Python`.
For F7, the raw `PA4` is also captured on CH2. For F4, the clock pin `PA8` (D7) is probed.

**Wiring notes**
- Op-amp supply: Li-ion cell straight to the supply pin, with 100 nF decoupling across the supply pins. The 470 Ω shunt resistor used in the power tests is **removed**.
- One common ground: MCU, battery negative, filter, probe ground clip. Keep the probe ground lead **short**. A long ground lead gave heavy ringing on the 8 MHz MCO signal and false threshold crossings.

---

## 4. Firmware (STM32 code)

All firmware is in [`firmware/`](firmware). The three files behave identically. They differ only in how the timer reload (ARR) is set and in whether the MCO pin is switched on.

| File | Use |
|---|---|
| [`firmware/F1_HSI.c`](firmware/F1_HSI.c) | Main test firmware. Preset for F1 (HSI). **ARR is computed at run time** from the real timer clock. Set the CONFIG block to run any test **except F4** |
| [`firmware/F4_test.c`](firmware/F4_test.c) | Same as above **plus MCO output on PA8** (`MCO_ENABLE 1`). Used for **F4 only** |
| [`firmware/freq_test.c`](firmware/freq_test.c) | Earlier version. Same waveform and clocks, but ARR is a fixed constant (79 for HSI, 719 for HSE). Header lists the tests it covers: F7, F1, F2, F3, F6 |

### 4.1 CONFIG block (the only lines you edit)

```c
#define CLOCK_SOURCE    CLOCK_HSI      /* CLOCK_HSI  or  CLOCK_HSE            */
#define TEST_FREQ_HZ    5000.0f        /* 1000 / 2000 / 5000 / 10000 only     */
#define KP              0.0f           /* 0.0f normally;  1.0f for F6 (PM)    */
#define MCO_ENABLE      1              /* F4_test.c only: 1 for F4, 0 otherwise */
```

| Test | CLOCK_SOURCE | TEST_FREQ_HZ | KP | Firmware file |
|---|---|---|---|---|
| F1 (T1) | HSI | 1000, 2000, 5000, 10000 | 0 | `F1_HSI.c` |
| F2 (T2) | HSE | 1000, 2000, 5000, 10000 | 0 | `F1_HSI.c` (CLOCK_HSE) or `freq_test.c` |
| F3 (T3) | HSI, then HSE (alternate) | 5000 | 0 | `F1_HSI.c` or `freq_test.c` |
| F4 (T4) | HSI, then HSE | 5000 (DAC runs but is not measured) | 0 | `F4_test.c` with `MCO_ENABLE 1` |
| F5 (T5) | HSI and HSE (two runs) | 5000 | 0 | `F1_HSI.c` or `freq_test.c` |
| F6 (T6) | HSE | 5000 | **1.0** | `F1_HSI.c` (KP 1.0f) or `freq_test.c` |
| F7 (T7) | HSE | 5000 | 0 | `F1_HSI.c` or `freq_test.c` |

### 4.2 What the firmware does

- Fills a 1000-sample buffer with a cosine, then streams it to the DAC by **TIM6 trigger → circular DMA**. The CPU does nothing during output.
- **Sample rate is exactly 100 kSa/s** on both clock paths.

| | HSI path | HSE path |
|---|---|---|
| Source | HSI 8 MHz **direct**, PLL off | HSE 8 MHz bypass (ST-LINK), PREDIV ÷1, PLL ×9 |
| SYSCLK | 8 MHz | 72 MHz |
| APB1 | ÷1 | ÷2 |
| TIM6 clock | 8 MHz | 72 MHz (timer doubling) |
| PSC / ARR | 0 / 79 | 0 / 719 |
| Flash latency | 0 | 2 |

- **Why the PLL is off on HSI.** On the F303xE the HAL constant `RCC_PLLSOURCE_HSI` can mean HSI/2 or HSI/PREDIV depending on the library version. With HSI/PREDIV a ×16 PLL gives 128 MHz, which is past the 72 MHz limit, and the DAC output collapses (this was observed). Running HSI directly avoids the problem. The core speed does not matter for this metric, because only the sample rate does.
- **Whole-cycle rule.** The buffer holds N = 1000 samples, so it repeats every 10 ms. Test frequencies must be multiples of 100 Hz so a whole number of carrier cycles fits (1, 2, 5, 10 kHz = 10, 20, 50, 100 cycles). Otherwise a phase jump appears at each DMA wrap.
- `DAC_MID = 2048`, `DAC_AMP = 1200` (about 0.68 V to 2.62 V).
- **Debug variables** (check in the debugger that the clock really switched): `dbg_sysclk`, `dbg_timclk`, `dbg_arr`, `dbg_ncycles` (and `dbg_mco` in `F4_test.c`). Expected values:

| | HSI | HSE |
|---|---|---|
| `dbg_sysclk` and `dbg_timclk` | 8 000 000 | 72 000 000 |
| `dbg_arr` | 79 | 719 |
| `dbg_ncycles` at 5 kHz | 50 | 50 |

> Do **not** regenerate the CubeMX code after pasting this `main.c`. It overwrites `SystemClock_Config()`. Keep the HAL timer in the `.ioc` as PSC = 0; the firmware sets ARR itself in `F1_HSI.c` and `F4_test.c`.

---

## 5. Oscilloscope settings

The scope always records **2000 points**, regardless of timebase or USB route. The sample interval therefore depends on the timebase, and that sets the quality of the estimate. The interval is in **row 2 of the CSV** (`Sample Interval`), and the Python code reads it from there.

Common settings: **AC coupling** (except F4), **Sample** acquisition, 20 MHz bandwidth limit ON, **Single** trigger on a rising edge at 0 V.

| Test | Probe point | Time/div | Sample interval | Volt/div | Notes |
|---|---|---|---|---|---|
| F1, F2, F3, F5, F6 at 5 kHz | LPF output, CH1 | **2 ms** | 16 µs (12.5 samples/cycle, about 160 cycles, about 1.1 ppm) | 500 mV | 200 mV/div sometimes clipped |
| F1, F2 at 1 kHz | LPF output | 10 ms | 80 µs | 500 mV | Frequency-matched timebase |
| F1, F2 at 2 kHz | LPF output | 5 ms | 40 µs | 500 mV | |
| F1, F2 at 10 kHz | LPF output | 1 ms | 8 µs | 500 mV | |
| F4 | **PA8 (MCO)** | **40 ns** | 1 GSa/s | 500 mV | **Bandwidth limit ON** (as in `t4.png`), **DC coupling**, trigger 0 V, shortest ground lead |
| F7 | CH1: LPF output, CH2: raw PA4 | 1 ms | 8 µs | 500 mV | One CSV holds both channels |

Do not use 10 ms/div at 5 kHz. It gives only 2.5 samples per cycle and a bias of about +5 ppm.

Scope images in [`oscilloscope_image/`](oscilloscope_image):

| File | Shows |
|---|---|
| `carrier_LPF.png` | Filtered 5 kHz carrier, about 1.88 Vpp |
| `raw_close.png` | Raw DAC staircase, close-up (20 µs/div) |
| `RAWvsFILTERED.png` | Filtered output (CH1) and raw DAC (CH2) together |
| `carrier_FFT.bmp` | Spectrum of the unmodulated carrier |
| `PM_FFT.bmp` | Spectrum with PM at β = 1 (carrier plus sidebands) |
| `t4.png` | MCO (PA8) 8 MHz clock output, 40 ns/div, 1 GSa/s, reads 8.032 MHz on the scope |

---

## 6. Folder layout and file-name mapping

```
frequency_accuracy/
├── firmware/            F1_HSI.c, F4_test.c, freq_test.c
├── python/              analyse_freq.py, batch_freq.py, plot_drift.py,
│                        analyse_stability.py, make_figures.py
├── csv/
│   ├── T1/              12 files + frequency_results.csv + "python result.txt"
│   ├── T2/              12 files + frequency_results.csv (this table also holds T3)
│   ├── T3/              6 files (analysed together with T2)
│   ├── T4/              6 files + frequency_results.csv
│   ├── T5_HSI/, T5_HSE/ 24 files each; T5_drift_results.csv, plots (in T5_HSE/)
│   ├── T6/              3 files + frequency_results.csv
│   └── T7/              T7.csv (CH1 + CH2 in one file)
├── oscilloscope_image/  scope screenshots
├── result_plots/        figA, figC, figD, figE, figF, T5_warmup_drift
└── Frequency_Accuracy_procedure.md   (this file)
```

### 6.1 File names were shortened in this repository

The scripts read the **test, clock and frequency from the file name**. The original names were as shown in the right-hand column. Several files in `csv/` were later renamed, so **rename them back before running the scripts** (this is done by the script in Section 7.1).

| Test | Name in `csv/` | Name the scripts expect |
|---|---|---|
| T1 (HSI) | `T1-1k1.csv`, `T1-1k2.csv`, `T1-1k3.csv` | `T1_HSI_1000Hz_run1.csv` … `run3.csv` |
| | `T1-2k1..3.csv`, `T1-5k1..3.csv`, `T1-10k1..3.csv` | `T1_HSI_2000Hz_run*`, `T1_HSI_5000Hz_run*`, `T1_HSI_10000Hz_run*` |
| T2 (HSE) | `T2_HSE_1000Hz_run1.csv` … | already correct |
| T3 | `T3_HSE_5000Hz_run*.csv`, `T3_HSI_5000Hz_run*.csv` | already correct |
| T4 | `HSE1.csv`, `HSE2.csv`, `HSE3.csv` | `T4_HSE_8MHz_run1.csv` … `run3.csv` |
| | `HSI1.csv`, `HSI2.csv`, `HSI3.csv` | `T4_HSI_8MHz_run1.csv` … `run3.csv` |
| T5 | `T5_HSI_t0030s.csv` … | already correct (the number is elapsed seconds) |
| T6 | `T6_HSE_5000Hz_run*.csv` | already correct |
| T7 | `T7.csv` | already correct (the frequency is passed on the command line) |

The renaming was checked: the renamed files give the same numbers as the saved `frequency_results.csv` files.

---

## 7. Reproducing the results

### 7.0 Setup

```bash
pip install "numpy>=2.0" matplotlib
```

(`analyse_stability.py` uses `numpy.trapezoid`, which needs numpy 2.0 or later.)

### 7.1 Stage the files (rename and group them)

Run this once from `frequency_accuracy/`. It copies everything into a new `work/` folder, so the repository files stay untouched. It also puts T2 and T3 together and both T5 clocks together, as in the original analysis.

```python
# stage.py  -- run from the frequency_accuracy folder:  python stage.py
import shutil, pathlib, re

repo = pathlib.Path(".")
work = pathlib.Path("work")
for d in ["T1", "T2_T3", "T4", "T5", "T6", "T7", "figs"]:
    (work / d).mkdir(parents=True, exist_ok=True)

for f in (repo / "python").glob("*.py"):                  # scripts into every folder
    for d in ["T1", "T2_T3", "T4", "T5", "T6", "T7", "figs"]:
        shutil.copy(f, work / d / f.name)

for f in (repo / "csv/T1").glob("T1-*.csv"):               # T1: restore long names
    m = re.match(r"T1-(\d+)k(\d)\.csv", f.name)
    shutil.copy(f, work / "T1" / f"T1_HSI_{int(m[1]) * 1000}Hz_run{m[2]}.csv")

for t, dst in [("T2", "T2_T3"), ("T3", "T2_T3"), ("T6", "T6")]:
    for f in (repo / "csv" / t).glob("T*_run*.csv"):
        shutil.copy(f, work / dst / f.name)

for f in (repo / "csv/T4").glob("HS*.csv"):                # T4: restore long names
    m = re.match(r"(HS[EI])(\d)\.csv", f.name)
    shutil.copy(f, work / "T4" / f"T4_{m[1]}_8MHz_run{m[2]}.csv")

for t in ["T5_HSI", "T5_HSE"]:
    for f in (repo / "csv" / t).glob("T5_*_t*.csv"):
        shutil.copy(f, work / "T5" / f.name)

shutil.copy(repo / "csv/T7/T7.csv", work / "T7" / "T7.csv")

for src in ["csv/T2/T2_HSE_5000Hz_run1.csv",               # inputs for make_figures.py
            "csv/T6/T6_HSE_5000Hz_run1.csv"]:
    shutil.copy(repo / src, work / "figs" / pathlib.Path(src).name)
shutil.copy(work / "T1/T1_HSI_5000Hz_run1.csv", work / "figs")
```

### 7.2 Analysis order and commands

The order below starts with the chain check (F7), because it justifies the filtered output as the measurement point.

| Step | Test | Command (run inside the folder shown) | Script | CSV files | Expected result |
|---|---|---|---|---|---|
| 1 | **F7 (T7)** chain check | `cd work/T7` then `python analyse_freq.py 5000 T7.csv` | `analyse_freq.py` | `T7.csv` (CH1 filtered, CH2 raw) | CH1 **+36.05 ppm**; CH2 +0.22 ppm, flagged as a staircase; verdict **PASS** (tolerance ±60 ppm) |
| 2 | **F1 (T1)** HSI sweep | `cd work/T1` then `python batch_freq.py` | `batch_freq.py` | 12 files: `T1_HSI_{1000,2000,5000,10000}Hz_run{1,2,3}.csv` | Mean **+3106.8 ppm**; per-frequency means 3190.5 / 3094.6 / 3086.9 / 3055.3 |
| 3 | **F4 (T4)** MCO pin | `cd work/T4` then `python batch_freq.py` | `batch_freq.py` | 6 files: `T4_{HSI,HSE}_8MHz_run{1,2,3}.csv` | HSI **+3047.2 ppm**; HSE +70.8 ppm |
| 4 | **F2 (T2)** and **F3 (T3)** HSE sweep and ratiometric | `cd work/T2_T3` then `python batch_freq.py` | `batch_freq.py` | 12 T2 files + 6 T3 files | T2 means 65.1 / 71.2 / 89.2 / 105.5 ppm; T3 HSI +3198.9, HSE +73.1, **difference +3125.7 ppm** |
| 5 | **F6 (T6)** PM mode | `cd work/T6` then `python batch_freq.py` | `batch_freq.py` (uses **Method C**, FFT, automatically for "T6" names) | 3 files `T6_HSE_5000Hz_run{1,2,3}.csv` | Mean **+72.06 ppm** (σ 29.3) |
| 6 | **F5 (T5)** warm-up and stability | `cd work/T5` then `python plot_drift.py` then `python analyse_stability.py` | `plot_drift.py`, then `analyse_stability.py` | 24 HSI + 24 HSE captures `T5_{HSI,HSE}_tNNNNs.csv` | HSI scatter 86.1 ppm, Allan 80.6, intrinsic 84.2; HSE 17.9 / 15.9 ppm. Trends not significant |
| 7 | Report figures | `cd work/figs` then `python make_figures.py` | `make_figures.py` | `T2_HSE_5000Hz_run1.csv`, `T6_HSE_5000Hz_run1.csv`, `T1_HSI_5000Hz_run1.csv` | `figures/figA, C, D, E, F` (PDF and PNG) |

Each `batch_freq.py` run prints a table, then writes `frequency_results.csv`. The copies in `csv/T1`, `csv/T4`, `csv/T6` and `csv/T2` are exactly those outputs (compared byte for byte). `csv/T2/frequency_results.csv` holds T2 and T3 together.

- `batch_freq.py` also prints a per-clock consistency check across frequencies and, when T1 and T4 data are in the same folder, an MCO-versus-DAC comparison. To get that comparison, copy the T1 files into `work/T4` and run `batch_freq.py` again.
- `plot_drift.py` writes `T5_drift_results.csv`, `T5_warmup_drift.png` and `.pdf`. `analyse_stability.py` reads `T5_drift_results.csv`, so run `plot_drift.py` first.

### 7.3 Figures

`make_figures.py` has a settings block at the top. **Figures A, D and F use numbers typed into that block** (the means and σ from `batch_freq.py` and `analyse_stability.py`), not numbers read from CSV files. Figures C and E read the CSV files directly. If you change the data, update the block.

| Figure | File | Shows |
|---|---|---|
| A | `figA_error_vs_frequency.png` | Error in ppm versus test frequency, HSI and HSE |
| C | `figC_spectrum_unmodulated_vs_pm.png` | Carrier spectrum, unmodulated versus PM (β = 1) |
| D | `figD_error_by_measurement_point.png` | HSI error at the MCO pin, raw DAC and filter output |
| E | `figE_cumulative_timing_error.png` | Cumulative timing error, HSI versus HSE |
| F | `figF_uncertainty_budget.png` | Uncertainty budget |
| – | `T5_warmup_drift.png` | HSI and HSE drift over about 31 min (from `plot_drift.py`) |

Regenerated PNGs match the committed plots in content, but may differ in the last bits of the file because of the matplotlib version.

---

## 8. Method details

All of this is in [`python/analyse_freq.py`](python/analyse_freq.py).

- **Loader (`load_tek_csv`).** Reads the Tektronix export. It takes the sample interval from the header (row 2) and handles one or two channels. It rejects unusable files (too few points, impossible span, error > 100 000 ppm).
- **Method D, zero-crossing regression** (all unmodulated tests). It finds the rising zero crossings with hysteresis (25 % of the amplitude) and interpolates between samples. It assigns cycle numbers from the crossing *time*, not by counting, so one missed crossing cannot shift every later index. It fits a straight line (crossing time against cycle number); the slope is the period. It trims 1 % of the record at each end and flags a bad fit.
- **Method C, interpolated FFT** (F6, PM mode). Modulation moves the zero crossings, so Method D is biased there (about 57 ppm at β = 1). Method C uses a Blackman-Harris window and parabolic interpolation of the spectral peak, which is not affected by this.
- **Staircase detection.** The raw DAC channel (a sample-and-hold staircase) is flagged. The crossings snap to the scope's sample grid (at 8 µs per scope sample and 10 µs per DAC step the two rates are commensurate, 5:4, which locks the fit to the nominal value), so its "σ" is far too optimistic and the true uncertainty is tens of ppm. **This is why the filtered output is the reported measurement point.** `batch_freq.py` shows such files as `STAIRCASE`.
- **PM check.** The β = 1 rad modulation index is confirmed from the sideband ratios: J1/J0 = −4.81 dB and J2/J0 = −16.47 dB for β = 1.

---

## 9. Known limitations and data notes

1. **Duplicate captures.** Four pairs of CSV files are byte-identical copies, probably from saving the same capture twice: T1 1 kHz run 1 and run 2 (+3256.83 ppm), T1 10 kHz run 2 and run 3 (+2999.76 ppm), T4 HSI run 1 and run 2 (+3107.89 ppm), and T2 10 kHz run 2 and run 3 (+136.12 ppm). They are kept in the data and counted as in the saved results. This slightly understates the true scatter for those groups. `plot_drift.py` checks for duplicates automatically (T5 has none).
2. **HSE is a bound, not a value.** The scope timebase (±20 ppm) is as large as the HSE tolerance. The F2 HSE results (about +65 to +105 ppm) include any scope timebase error. The scope also shows a trend with timebase range (about +4.5 ppm per kHz), because different frequencies use different time/div settings. The 10 kHz HSE group is the noisiest (σ 53 ppm, with a duplicate pair).
3. **F5 timing.** The first capture is at **30 s (HSI)** and **60 s (HSE)** after power-up, not at t = 0. Capture gaps are 30 s at first and then 2 min. So "drift at t = 0" in the output is an extrapolation and not a measurement. No significant warm-up trend was found.
4. **2000-point record.** The scope's record length is fixed. This limits how many cycles can be captured and sets the estimator precision (about 1.1 ppm at 2 ms/div, 5 kHz).
5. **Probe setting.** Probe and scope channel must both be set to 10×.
6. **Test names in the code.** `batch_freq.py` finds the test number with the pattern `T1`–`T7` in the file name. It does not recognise `F1`–`F7`. If you rename the CSVs to F-names, change that pattern in `test_and_clock()`.
7. **Firmware header comments.** `F4_test.c` still says the MCO should be measured with the bandwidth limit OFF. The measurement in `t4.png` was taken with it ON. Section 5 lists the setting used for the recorded data.

---

## 10. Quick start (no hardware)

```bash
git clone https://github.com/moyeenurrahman/STM32_phase_modulation_research
cd STM32_phase_modulation_research/frequency_accuracy
pip install "numpy>=2.0" matplotlib
python stage.py          # save the script from Section 7.1 as stage.py first
cd work/T1 && python batch_freq.py     # HSI sweep, expect mean about +3106.8 ppm
```

Then follow the table in Section 7.2 for the other tests.
