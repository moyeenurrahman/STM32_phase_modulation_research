# SNR, SINAD and ENOB Measurement: Procedure and Reproduction Guide

Part of the thesis *"Implementation and Performance Analysis of a Baseband Phase Modulation (PM) System with STM32 for Low-Power Wireless Communication"* (BSc EEE, RUET).

This guide explains how the **signal-to-noise ratio (SNR)**, **SINAD** and **ENOB** of the STM32 baseband signal generator were measured. It tells you which firmware was flashed, how the oscilloscope was set, which CSV files were recorded, which Python script turns each set of CSV files into a result, and what each plot and oscilloscope image shows.

You can follow it to **reproduce every number and plot in this folder from the CSV files alone, without any hardware.** Section 3 does this in about five minutes.

---

## 1. What is measured, and the main results

The STM32 makes a 5 kHz carrier with its 12-bit DAC. The signal is measured at two places:

- **Test Point A:** the raw DAC output (a staircase).
- **Test Point B:** the output of a Sallen-Key low-pass filter (the clean signal).

Six cases (N1 to N6) were captured. Two analysis methods were used:

| Method | Script | What it gives |
|---|---|---|
| **A: spectrum method** | `snr_method_a.py` | In-band SNR (0.1 to 50 kHz), wideband SNR, SFDR, and the DAC image levels |
| **B: sine-fit method** (IEEE Std 1241) | `snr_method_b_sinad.py` | SINAD and ENOB, plus a fitted carrier frequency |

### Headline results (all reproducible from this folder)

| Case | Signal | Test point | In-band SNR, standard | In-band SNR, robust | Wideband SNR | SFDR |
|---|---|---|---|---|---|---|
| N1 | Noise floor (DAC at mid-code) | B | −14.5 dB (floor, not a result) | n/a | n/a | n/a |
| N2 | 5 kHz carrier | A (raw DAC) | **29.11 dB** | 39.82 dB | 20.76 dB | 25.59 dB |
| N3 | 5 kHz carrier | B (filtered) | **38.62 dB** | 41.09 dB | 32.40 dB | 42.21 dB |
| N4 | PM, kp = 1.0 | B | 32.18 dB (band 3.5 to 6.5 kHz) | 38.82 dB | 29.94 dB | n/a |
| N5 | PM, kp = 2.4048 | B | 30.45 dB (band 2.8 to 7.2 kHz) | 38.30 dB | 28.65 dB | n/a |
| N6 | BPSK, 500 bps | B | 18.71 dB (band 1.25 to 8.75 kHz) | 35.86 dB | 19.24 dB | n/a |

| Quantity | Point A (N2) | Point B (N3) | Change across the filter |
|---|---|---|---|
| SINAD (Method B) | 20.69 dB | 32.32 dB | +11.63 dB |
| ENOB (Method B) | 3.92 bits | 5.85 bits | +1.93 bits |
| Wideband SNR (Method A) | 20.76 dB | 32.40 dB | +11.64 dB (agrees with SINAD) |
| Image at 95 kHz | −25.59 dBc | −42.21 dBc | 16.6 dB lower |
| Image at 105 kHz | −26.56 dBc | −43.47 dBc | 16.9 dB lower |
| Fitted carrier frequency | 5016.091 Hz (σ 0.350 Hz) | 5016.187 Hz (σ 0.520 Hz) | n/a |

### How to read these results honestly

1. **The measurement is limited by the oscilloscope and by interference, not by the DAC.** An ideal 12-bit DAC at this amplitude would give about 69 dB (Point A) and 73 dB (Point B). The 8-bit oscilloscope can show at most about 44.8 dB in this setup (the "scope limit"). The measured values are below that limit because discrete interference lines add to the noise.
2. **Standard and robust SNR differ because of interference lines.** The *standard* value counts everything in the band as noise. The *robust* value drops bins that are more than 10 times the median (single spurs, such as pick-up from other circuits). A big gap between the two means a few discrete lines dominate. For N2, 28 bins were dropped. For N3, 11.
3. **N4 to N6 are not IEEE 1241 SNR.** The sidebands of a phase-modulated signal have no sharp edge, so the signal band is a choice. These values are "in-band SNR with the band stated", and are useful only for comparing conditions with each other.
4. **ENOB here describes the whole measurement chain** (DAC, filter, scope and interference together). It is not the ENOB of the DAC alone.
5. The fitted carrier frequency (5016 Hz, about +0.32 % high) is consistent with the STM32 **HSI** internal oscillator error found in the separate `frequency_accuracy` folder (+3107 ppm there, +3218 ppm here). The two are independent measurements on different days, and the HSI drifts with temperature and supply, so exact agreement is not expected.

---

## 2. Folder map

```
SNR/
├── SNR_procedure.md              <- this guide
├── firmware/
│   └── SNR.c                     <- the only firmware file (copy it into a CubeIDE project as main.c)
├── csv/                          <- oscilloscope records, 2000 points each
│   ├── N1/   32 files (1.csv ... 32.csv)  + "python result.txt"
│   ├── N2/   32 files                     + "python result method A.txt", "python result method B.txt"
│   ├── N3/   32 files                     + "python result method A.txt", "python result method B.txt"
│   ├── N3_configB/  16 files
│   ├── N4/   8 files + "python result.txt"
│   ├── N5/   8 files + "python result.txt"
│   └── N6/   8 files + "python result.txt"
├── python/
│   ├── check_csv.py              <- checks ONE csv file (sample rate, clipping, aliasing)
│   ├── snr_method_a.py           <- Method A (spectrum): SNR, SFDR, images
│   ├── snr_method_b_sinad.py     <- Method B (sine fit): SINAD, ENOB, frequency
│   ├── measure_b99.py            <- 99 % power bandwidth of the BPSK signal (N6)
│   └── make_summary_figures.py   <- three summary bar charts (no CSV needed)
├── result_plots/                 <- plots made by the scripts (Section 8)
└── oscilloscope_images/          <- scope screenshots (Section 9)
```

Every CSV folder also holds the text output the script printed when it was run. You can compare your own run against these.

---

## 3. Quick reproduction (no hardware)

**Requirements:** Python 3.9 or newer with `numpy`, `scipy` and `matplotlib`.

```bash
pip install numpy scipy matplotlib
```

All scripts have a **configuration block at the top**. You edit it, then run the script. The plots are saved **in the folder you run the script from**.

Example: reproduce case N3.

1. Open a terminal in `SNR/python/`.
2. Open `snr_method_a.py` and set:
   ```python
   FOLDER       = "../csv/N3"
   FLOOR_FOLDER = "../csv/N1"
   LABEL        = "N3 - filtered output"
   SIGNAL_MODE  = "tone"
   ```
3. Run `python snr_method_a.py`.
4. You should see `In-band SNR (0.1-50 kHz) : 38.62 dB`, `robust variant : 41.09 dB` and `Wideband SNR (0.1-148 kHz): 32.40 dB`. The plots `psd_N3.png` and `images_N3.png` are written to the folder you ran from.

Section 7 lists the exact settings for every run.

---

## 4. The test cases

| Case | What it is | Firmware mode (`SNR_MODE`) | Probe at | CSV folder | Records | Script |
|---|---|---|---|---|---|---|
| **N1** | Instrument noise floor. The DAC output is constant (mid-code), so any signal seen is scope and pick-up noise | `MODE_DC_DMA` (0) | Test Point B | `csv/N1` | 32 | A (floor for all others) |
| **N2** | 5 kHz carrier, raw DAC output | `MODE_CARRIER` (1) | Test Point A | `csv/N2` | 32 | A and B |
| **N3** | 5 kHz carrier, filtered output | `MODE_CARRIER` (1) | Test Point B | `csv/N3` | 32 | A and B |
| **N3 config B** | N3 repeated with a different power arrangement | `MODE_CARRIER` (1) | Test Point B | `csv/N3_configB` | 16 | A |
| **N4** | Analog PM, phase deviation kp = 1.0, fm = 500 Hz | `MODE_PM_1` (2) | Test Point B | `csv/N4` | 8 | A (pm mode) |
| **N5** | Analog PM, kp = 2.4048 (carrier null of the Bessel function J0) | `MODE_PM_BESSEL` (3) | Test Point B | `csv/N5` | 8 | A (pm mode) |
| **N6** | BPSK, 500 bps, 5-bit repeating pattern | `MODE_BPSK` (4) | Test Point B | `csv/N6` | 8 | `measure_b99.py`, then A (bpsk mode) |

**Power arrangements for N3.** *Config A* (folder `N3`): STM32 powered from a phone charger, laptop on battery. *Config B* (folder `N3_configB`): STM32 powered from the laptop USB, laptop on battery. All captures were taken with the laptop on its battery. A third arrangement (STM32 on laptop USB with the laptop on its mains charger) showed visible oscillation on the scope, so no data set was recorded for it.

**Why a noise-floor case (N1).** Without it you cannot tell whether a measured SNR belongs to the signal or to the oscilloscope. Every analysis run compares against `csv/N1`, and the script prints whether the result is a *measured value* or only a *lower bound*.

---

## 5. Hardware and firmware

### 5.1 Signal chain

```
STM32 Nucleo-F303RE
  TIM6 TRGO --> DAC1_OUT1 (PA4), circular DMA, 1000-sample buffer
                    |
                    +--> Test Point A  (raw DAC staircase)
                    |
              Sallen-Key 2nd-order Butterworth low-pass filter
              R1 = R2 = 5.6 kOhm, C1 = 2 nF, C2 = 1 nF (fc about 20.1 kHz), MCP6002 op-amp
                    |
                    +--> Test Point B  (filtered output)
```

### 5.2 Fixed parameters

| Parameter | Value |
|---|---|
| DAC sample rate fs | 100 kSa/s (TIM6 on an 8 MHz clock, ARR = 79) |
| Carrier fc | 5 kHz (= fs / 20) |
| Modulating frequency fm (PM modes) | 500 Hz |
| BPSK symbol rate | 500 bps (200 DAC samples per symbol) |
| DMA buffer | 1000 samples = exactly 50 carrier cycles, 5 modulating cycles, 5 BPSK symbols |
| DAC mid-level, amplitude | 2048 codes, 1200 codes (about 0.68 to 2.62 V, about 1.94 Vpp) |
| Clock | HSI at 8 MHz, no PLL (`SystemClock_Config` in `SNR.c`) |

The buffer holds a **whole number** of carrier, modulation and symbol cycles. If it did not, the waveform would jump every time the DMA wraps around (every 10 ms), and this would add a comb of lines 100 Hz apart to the spectrum. The firmware checks this at compile time and refuses to build otherwise.

### 5.3 The firmware: one file, six builds

`firmware/SNR.c` is a complete `main.c` for an STM32CubeIDE project (Nucleo-F303RE, with DAC1, DMA, TIM6 and GPIO generated by CubeMX). **One line selects the build:**

```c
#define SNR_MODE   MODE_CARRIER
```

| `SNR_MODE` | LED flashes at start | Output | Used for |
|---|---|---|---|
| `MODE_DC_DMA` (0) | 1 | DAC at mid-code (2048), DMA running | N1 |
| `MODE_CARRIER` (1) | 2 | Pure 5 kHz sine | N2, N3 |
| `MODE_PM_1` (2) | 3 | Analog PM, kp = 1.0 | N4 |
| `MODE_PM_BESSEL` (3) | 4 | Analog PM, kp = 2.4048 | N5 |
| `MODE_BPSK` (4) | 5 | BPSK 500 bps, 5-bit pattern `0x16` (binary 10110) repeating every buffer | N6 |
| `MODE_DC_STATIC` (5) | 6 | DMA stopped, DAC held at 2048 | Optional. Not used in the final data |

The LED flash count at power-up tells you which build is loaded, without a debugger.

**Design choices that matter for a noise measurement:**

- The main `while(1)` loop is **empty**. No printing, no LED toggling, no sleep mode. Anything the CPU does would add switching noise to the supply and appear in the spectrum.
- `force_dac_config()`, `force_dma_config()` and `force_sample_rate()` re-set the settings the measurement depends on, so a regenerated CubeMX project cannot change them silently.
- `MCO_ENABLE` must be 0. The clock output used in the frequency accuracy test would put a tone on the supply rails. The build stops with an error if it is on.
- `Error_Handler()` makes the LED blink fast, so a failed start is visible.

### 5.4 Flashing

1. Create a CubeIDE project for the Nucleo-F303RE with DAC1 (OUT1 on PA4, DMA circular), TIM6 (TRGO update event) and GPIO. Keep CubeMX **peripheral init as separate files** (`dac.c`, `dma.c`, `gpio.c`, `tim.c`).
2. Copy `firmware/SNR.c` over the project's `Core/Src/main.c`.
3. Set `SNR_MODE`, then **Project → Clean**, **Project → Build**.
4. Flash with **Run** (not Debug).
5. **Unplug the USB** and power the board from a charger (or USB power), as in Section 6.2. The LED flashes the mode number once to confirm the build.

---

## 6. Capture procedure (oscilloscope)

### 6.1 Equipment

Tektronix **TBS1102C** oscilloscope (8-bit, 2000-point records), 1× passive probe, laptop with **TekVISA** and **OpenChoice Desktop**, USB cable (scope USB device port to laptop), the STM32 board and the filter board.

### 6.2 Wiring and power

1. Keep wires between the Nucleo and the filter board short (under 10 cm).
2. PA4 goes to the filter input (this node is **Test Point A**). The filter output is **Test Point B**.
3. Use the **shortest probe ground lead** on the shared ground. A long lead picks up noise, and in an SNR test the pick-up becomes part of the result.
4. Run the laptop **on battery** during capture. With the laptop on its mains charger, visible oscillation appeared on the scope.
5. Move phone chargers, soldering irons and other switching supplies away from the bench.

### 6.3 Oscilloscope settings (the same for every case)

| Setting | Value | Reason |
|---|---|---|
| Channel | CH1 | n/a |
| Probe | 1× (switch **and** menu) | Full amplitude to the 8-bit ADC |
| Coupling | **AC** | Removes the DC offset so the signal fills the screen |
| Vertical scale | **500 mV/div** | The 1.94 Vpp signal fits without clipping |
| Acquisition | **Sample** | **Not** Hi-Res (see below) |
| Averaging | **Off** | Averaging would remove the noise you are measuring and give a falsely high SNR |
| Bandwidth limit | 20 MHz on | Cuts front-end noise above 20 MHz |
| Timebase | **400 µs/div** | 312.5 kSa/s, 2000 points, 6.4 ms per record |
| Trigger | Edge, CH1, 0 V, **Normal** (use **Auto** for N1, which has no signal) | n/a |

The vertical scale must **not change between Point A and Point B**. Only the probe tip is moved, so that both can be compared on equal terms.

What 400 µs/div gives:

| Quantity | Value |
|---|---|
| Sample rate | 312.5 kSa/s (rule: fs = 125 / timebase in ms, in kSa/s) |
| Record length | 2000 points = 6.4 ms |
| FFT bin width | 156.25 Hz (= 62.5 / timebase in ms, in Hz) |
| Nyquist | 156.25 kHz |

**Do not use 1 ms/div.** At 1 ms/div the scope samples at 125 kSa/s, and the 95 kHz DAC image would fold to 30 kHz, inside the noise band. The scripts also stop with an error if an image folds into the 0.1 to 50 kHz band.

**Do not use Hi-Res mode.** Hi-Res returns 16-bit values, and OpenChoice then saves a wrong waveform. Sample mode saves the true 8-bit data.

### 6.4 Saving records with OpenChoice

1. Set the scope USB port to **Computer** (Utility → Options → USB Device Port).
2. In OpenChoice Desktop, select the instrument.
3. Press **Single** on the scope and **wait until it shows Stop**.
4. In OpenChoice: **Waveform Data Capture → CH1 → Get Data → Save as CSV**.
5. Repeat for 32 records (N1, N2, N3), 16 (N3 config B) or 8 (N4, N5, N6). Name the files `1.csv`, `2.csv`, and so on, one folder per case.

If you click **Get Data** before the scope has stopped, you get a smeared record (a "solid band"). Always wait for Stop.

### 6.5 Order of work

1. Flash `MODE_DC_DMA`, capture **N1** (trigger on Auto).
2. Flash `MODE_CARRIER`. Capture **N2** (probe at Point A), then **N3** (probe at Point B).
3. Capture **N3 config B** with the other power arrangement.
4. Flash `MODE_PM_1`, `MODE_PM_BESSEL`, `MODE_BPSK` in turn and capture **N4, N5, N6** at Point B.
5. Take scope screenshots (Section 9) for each case.

### 6.6 CSV file format

Each file is one 2000-point record in the Tektronix format. The first rows carry the settings in the first columns, and the data is in the last two columns:

```
Record Length,2000,Points,-3.2E-3,0.0,
Sample Interval,3.2E-6,s,-3.1968E-3,0.0,
...
Vertical Scale,500E-3,...
```

The last two columns are time (s) and voltage (V). The scripts read the sample rate from the *Sample Interval* line (3.2 µs gives 312.5 kSa/s). Files with only a voltage column are also accepted, using `FS_FALLBACK` (default 312.5 kSa/s).

---

## 7. Analysis: which script, which CSV folder, which settings

Run the scripts from `SNR/python/`. The plots are saved in the folder you run from. Edit only the lines shown.

### 7.1 Order of runs

1. (Optional) `check_csv.py` on one file, to confirm the record is healthy.
2. `snr_method_a.py` on **N1** against itself. This confirms the floor.
3. `snr_method_a.py` on N2, N3, N3 config B.
4. `snr_method_b_sinad.py` on N2 and N3.
5. `snr_method_a.py` on N4 and N5 (pm mode).
6. `measure_b99.py` on N6, **then** `snr_method_a.py` on N6 (bpsk mode, using the B99% result).
7. `make_summary_figures.py` for the three summary charts.

### 7.2 Method A: `snr_method_a.py`

Settings to edit (all other settings keep their defaults):

| Run | `FOLDER` | `FLOOR_FOLDER` | `LABEL` | `SIGNAL_MODE` | `KP` | `BPSK_HALFBW` |
|---|---|---|---|---|---|---|
| N1 | `"../csv/N1"` | `"../csv/N1"` | `"N1 - instrument noise floor"` | `"tone"` | n/a | n/a |
| N2 | `"../csv/N2"` | `"../csv/N1"` | `"N2 - Raw DAC output"` | `"tone"` | n/a | n/a |
| N3 | `"../csv/N3"` | `"../csv/N1"` | `"N3 - filtered output"` | `"tone"` | n/a | n/a |
| N3 config B | `"../csv/N3_configB"` | `"../csv/N1"` | `"N3 - config B"` | `"tone"` | n/a | n/a |
| N4 | `"../csv/N4"` | `"../csv/N1"` | `"N4 - PM, kp=1.0"` | `"pm"` | `1.0` | n/a |
| N5 | `"../csv/N5"` | `"../csv/N1"` | `"N5 - PM, kp=2.4048"` | `"pm"` | `2.4048` | n/a |
| N6 | `"../csv/N6"` | `"../csv/N1"` | `"N6 - BPSK 500 bps"` | `"bpsk"` | n/a | `3750.0` |

Other settings: `FC = 5000`, `FM = 500`, `FS_FALLBACK = 312.5e3`. The wideband limit (`WIDE_TOP`) reduces itself automatically to 148 kHz, because the Nyquist limit is 156 kHz.

**What the script does:**

1. Reads every CSV in the folder (each file is read once, even if the folder lists it under two spellings).
2. Takes a power spectrum of each 2000-point record with a Blackman-Harris window, then **averages the spectra of all records**. Averaging 32 records of 2000 points reduces the scatter of the noise estimate in the same way as averaging 32 Welch segments in a single long record. (The oscilloscope can store only 2000 points, so one long record is not possible.)
3. Marks bins for the **signal** (±5 bins around the carrier), the **harmonics** 2 to 10 (they alias onto multiples of 5 kHz, because fc = fs/20), the **DAC images** (95 and 105 kHz) and **DC** (below 100 Hz). These are left out of the noise.
4. Noise power = average noise density × band width, for 0.1 to 50 kHz (in-band) and 0.1 to 148 kHz (wideband).
5. Gives a **standard** SNR and a **robust** SNR (bins above 10 times the median are dropped).
6. Compares the result with the N1 floor and prints a verdict: *MEASURED VALUE* if the noise is at least 3 dB above the floor, otherwise *LOWER BOUND*.

**Mode details:**

- `"tone"`: the IEEE 1241 single-tone SNR.
- `"pm"`: the signal band is (kp + 2) × fm each side of the carrier (N4: 3.5 to 6.5 kHz, N5: 2.8 to 7.2 kHz), wider than the strict Carson band. Harmonics are not excluded, because they cannot be told apart from modulation products.
- `"bpsk"`: the signal band is ±`BPSK_HALFBW` around the carrier. Use the value from `measure_b99.py`.

**Expected output:** the files `csv/N*/python result*.txt` hold the printed output of the original runs. Section 10 lists the key values.

### 7.3 Method B: `snr_method_b_sinad.py`

| Run | `FOLDER` | `LABEL` |
|---|---|---|
| N2 | `"../csv/N2"` | `"N2 - raw DAC output"` |
| N3 | `"../csv/N3"` | `"N3 - filtered output"` |

Other settings: `FS_DAC = 100000`, `DAC_CODES = 2400` (2 × 1200 codes used), `DAC_FULL = 4096`.

**What it does:** for each record it fits a sine (amplitude, phase, offset and **frequency**, i.e. a four-parameter fit, IEEE Std 1241). The leftover signal is the residual. SINAD = sine power / residual power. ENOB comes from SINAD, corrected for the fact that only 2400 of the 4096 DAC codes are used (4.64 dB):

`ENOB = (SINAD − 1.76 + 4.64) / 6.02`

It prints the mean and standard deviation across the 32 records, the fitted frequency, and a check that the residual looks like noise.

**Why SINAD is lower than the in-band SNR.** SINAD counts harmonic distortion and noise over the whole 0 to 156 kHz capture, while the in-band SNR counts only 0.1 to 50 kHz. If SINAD were higher than the in-band SNR, one of the analyses would be wrong.

### 7.4 B99%: `measure_b99.py` (BPSK only)

| Setting | Value |
|---|---|
| `FOLDER` | `"../csv/N6"` |
| `FC` | `5000.0` |
| `SEARCH_HALFBW` | `40000.0` (default) |

It grows a band around the carrier until the band holds 99 % of the signal power. It also repeats the test with a doubled window to show the answer has converged.

**Result:** B99% = **7.500 kHz** (±3.750 kHz). 99.9 % needs 12.8 kHz. This is why `BPSK_HALFBW = 3750.0` is used for N6 in Method A.

### 7.5 Checking one file: `check_csv.py`

Set `FILENAME = "1"` and place the script in the same folder as the CSV (or give a full path). It prints the sample rate, checks that the record is coherent (a clean carrier, not a smeared capture), finds the carrier frequency, looks for clipping and aliasing, and saves `check_csv.png`. Use it right after the first capture of a session to catch a wrong setting early.

### 7.6 Summary charts: `make_summary_figures.py`

This script does **not** read CSV files. It plots numbers typed into its `RESULTS` block (copied from the run outputs above), and writes `fig4_inband_snr_summary.png`, `fig6_sinad_enob_improvement.png` and `fig10_power_config_comparison.png`. If you re-run the analysis, check that these numbers still match.

---

## 8. Result plots (`result_plots/`)

All plots are made by the scripts in Section 7. The file name tells you the case.

### 8.1 Spectrum plots: `snr_method_a.py`

| File | What it shows | What to look for |
|---|---|---|
| `psd_N1.png` | Spectrum (power density, dB/Hz) of the **noise floor** case, 0 to 148 kHz | A flat line at about −95 dB/Hz. This is the oscilloscope floor. |
| `psd_N2.png` | Spectrum of the **raw DAC output** with the N1 floor overlaid. The green band is the signal band. The blue-grey region is the in-band noise region. Red dashed lines mark 95 and 105 kHz. | A sharp 5 kHz peak. Many narrow spurs well above the floor. Two large peaks at 95 and 105 kHz (DAC images, about −55 dB/Hz). |
| `psd_N3.png` | Same for the **filtered output** | The 5 kHz peak is the same height, the images are about 17 dB lower, and the spurs above 20 kHz are lower. |
| `psd_N4.png` | PM, kp = 1.0, with the signal band (3.5 to 6.5 kHz) | The carrier with sidebands around it. |
| `psd_N5.png` | PM, kp = 2.4048, signal band 2.8 to 7.2 kHz | Carrier strongly suppressed (−37.9 dB below the strongest line, near the Bessel null). |
| `psd_N6.png` | BPSK, signal band 1.25 to 8.75 kHz | A wide, smooth spectrum with a slow fall-off. **The title says "N3 - congigB" and "8 records". This is a leftover label. The data is the N6 data (8 records).** |
| `images_N1.png` to `images_N6.png` | Zoom on 80 to 120 kHz, the **DAC images** at 95 and 105 kHz, with the N1 floor | In `images_N2.png` the images stand about 40 dB above the floor. In `images_N3.png` they are about 17 dB lower. |

### 8.2 Sine-fit plots: `snr_method_b_sinad.py`

| File | What it shows | What to look for |
|---|---|---|
| `sinefit_N2.png`, `sinefit_N3.png` | Top: one measured record (blue) and the fitted sine (orange). Bottom: the residual (measured minus fit). | The residual should look like noise, not like a sine. In `sinefit_N2.png` the residual is larger (the staircase and interference). In `sinefit_N3.png` it is smaller. **`sinefit_N2.png` carries the title "N3 - raw DAC output". This is a leftover label; the data is N2.** |
| `residual_spectrum_N2.png`, `residual_spectrum_N3.png` | Spectrum of the residual | No special peak at the carrier (the script reports 3.9 dB above the median for N2: acceptable). The large peaks at 95 and 105 kHz and above are the DAC images and interference, which the fit leaves in the residual. |

### 8.3 BPSK bandwidth: `measure_b99.py`

| File | What it shows |
|---|---|
| `b99_N6.png` | Top: BPSK spectrum with the 99 % band shaded (1.25 to 8.75 kHz). Bottom: cumulative power against distance from the carrier. The curve crosses 99 % at ±3.75 kHz. |

### 8.4 Summary charts: `make_summary_figures.py`

| File | What it shows |
|---|---|
| `fig4_inband_snr_summary.png` | Standard and robust in-band SNR for N2 to N6, with the scope limit (44.8 dB) as a dashed line. N4 to N6 use a band that depends on the modulation, so they are not on the same footing as N2 and N3. |
| `fig6_sinad_enob_improvement.png` | SINAD (20.69 to 32.32 dB) and ENOB (3.92 to 5.85 bits) at Point A and Point B. The note under the chart explains that a gain this large cannot come from DAC quantization noise, which the filter does not remove below 20 kHz. It comes from interference above the filter cut-off. |
| `fig10_power_config_comparison.png` | N3 for power arrangement A (charger) and B (laptop USB): SNR 38.62 / 41.09 dB against 38.57 / 41.21 dB, and 11 against 9 dropped bins. The two arrangements are not different in a meaningful way. |

---

## 9. Oscilloscope images (`oscilloscope_images/`)

Screenshots from the TBS1102C. They are for illustration; the analysis uses the CSV files.

| File | What it shows |
|---|---|
| `N2N3.png` | Both channels at once, 40 µs/div. **CH1 (yellow): filtered sine** (Point B). **CH2 (cyan): raw DAC staircase** (Point A). It shows what the filter does in the time domain. |
| `N2fft.png` | CH2 (Point A) with the scope's own FFT view (20 dB/div, 400 µs/div). Many spurs across the span, plus the carrier at the left edge. |
| `N3.png` | CH1 (Point B) at 400 µs/div: a clean 5 kHz sine, 1.96 Vpp. This is the setting used for the CSV records. |
| `N3fft.png` | CH1 (Point B) FFT: a single carrier peak, a lower and cleaner floor than `N2fft.png`. |
| `N4_scope.png` | PM, kp = 1.0, at 100 µs/div. Overlapping traces: the phase changes from one sweep to the next, so a 0 V trigger shows several phases at once. This is the visible effect of phase modulation. |
| `N5_scope.png` | PM, kp = 2.4048, at 400 µs/div. The amplitude stays constant while the spacing of the cycles changes strongly: a large phase swing. |
| `N6_scope.png` | BPSK at 400 µs/div. The abrupt **phase reversals** at symbol boundaries show as distorted cycles in the carrier. |

The trigger-frequency readout in the bottom right of `N2N3.png`, `N2fft.png`, `N3.png` and `N3fft.png` shows 5.016 kHz, which agrees with the 5016 Hz fitted in Method B. The on-screen measurement boxes (for example 5.027 kHz) are single-shot readings and are less accurate.

---

## 10. Reference values for checking your run

| Run | Output to check |
|---|---|
| N1 | Floor noise power 1.5422e-05 V², about −14.5 dB. Verdict: LOWER BOUND, instrument limited |
| N2 (A) | In-band 29.11 dB, robust 39.82, wideband 20.76, SFDR 25.59 (at 95.31 kHz), images −25.59 / −26.56 dBc, 28 outlier bins |
| N3 (A) | In-band 38.62 dB, robust 41.09, wideband 32.40, SFDR 42.21, images −42.21 / −43.47 dBc, 11 outlier bins |
| N3 config B | In-band 38.57 dB, robust 41.21, 9 outlier bins |
| N4 | Band 3500 to 6500 Hz, in-band 32.18 dB, robust 38.82, wideband 29.94 |
| N5 | Band 2798 to 7202 Hz, in-band 30.45 dB, robust 38.30, wideband 28.65; strongest line 4531.2 Hz, carrier −37.9 dB |
| N6 | B99% = 7.500 kHz; band 1250 to 8750 Hz, in-band 18.71 dB, robust 35.86, wideband 19.24, 63 outlier bins |
| Method B, N2 | f0 5016.091 Hz, SINAD 20.69 dB (σ 0.03), ENOB 3.92 bits |
| Method B, N3 | f0 5016.187 Hz, SINAD 32.32 dB (σ 0.12), ENOB 5.85 bits |

All of the values above were reproduced by re-running the scripts in this folder on the CSV files in this folder (Python 3 with numpy 2.5.3, scipy 1.18.1, matplotlib 3.11.2). Small differences in the last digit between a re-run and a saved text file are rounding.

---

## 11. Known limits and notes

1. **Record counts in older text files.** In `csv/N1`, `csv/N2` and `csv/N3`, the saved "python result" files say *Records averaged: 64*. Each folder holds 32 files. An earlier version of the script read every file twice on Windows (once for `*.csv` and once for `*.CSV`). Averaging each record twice gives the **same** mean spectrum, so the SNR numbers are unchanged, but the record count printed was wrong. The current scripts count each file once and print 32. The N4 to N6 text files already show the correct floor count (32).
2. **Leftover plot titles.** `psd_N6.png` shows "N3 - congigB" and `sinefit_N2.png` shows "N3 - raw DAC output". Both labels are wrong; the data in each plot is correct (N6 and N2). The `python result.txt` in `csv/N6` has the same stale label in its header.
3. **Instrument limit.** An 8-bit oscilloscope cannot resolve the 69 to 73 dB predicted for an ideal 12-bit DAC. In-band SNR is therefore a measurement of the whole chain, with a ceiling of about 44.8 dB.
4. **Interference.** Most of the in-band noise at Point A sits in a few dozen discrete bins (not random noise). The low-pass filter removes much of it, which is why SINAD and wideband SNR improve by about 11.6 dB across the filter. The source of this interference was not isolated (possible causes include the ST-LINK supply and the lab environment).
5. **Image suppression.** The measured drop in image level across the filter is about 16.6 dB. The filter model alone would predict a larger drop, so the filter response and the interference floor probably limit it. The filter's frequency response was **not measured** in this work.
6. **Power arrangement C** (STM32 on laptop USB with the laptop on its mains charger) showed oscillation on the scope and was **not** quantified.
7. **N4 to N6** are not IEEE 1241 SNR. SFDR is not meaningful for them, because their spectra have no single dominant line.
8. **Summary charts** (`make_summary_figures.py`) use hand-entered numbers. They match the run outputs above.
9. **The firmware clock** is the internal HSI at 8 MHz (no PLL), the same oscillator as in the frequency accuracy measurement. This is why the carrier appears at 5016 Hz and not 5000 Hz.

---

## 12. Related folders

| Folder | Link to the SNR work |
|---|---|
| `frequency_accuracy/` | The HSI frequency error (+0.31 %) explains the 5016 Hz carrier. |
| `THD/` | Uses the same signal chain and a related analysis. |
| `phase_accuracy/` | The N5 carrier null (−37.9 dB) can be compared with the phase-chapter result near kp = 2.4073. |
| `power_consumption/` | The power measurements that raised the op-amp supply question (Section 11, item 6). |
