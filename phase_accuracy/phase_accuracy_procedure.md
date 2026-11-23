# Phase Deviation and Phase Accuracy: Measurement and Reproduction Guide

Part of the thesis *"Implementation and Performance Analysis of a Baseband Phase Modulation (PM) System with STM32 for Low-Power Wireless Communication"* (BSc EEE, RUET).

This guide explains how the **phase deviation and phase accuracy** of the STM32 baseband PM generator were measured. It tells you which firmware was flashed, how the oscilloscope was set, which CSV files were recorded, and which Python script turns each set of CSV files into a result. You can follow it to **reproduce every number and plot from the CSV files in this folder**, without any hardware.

> **Test names.** The tests were first called **M1, M2, M3** and later renamed **PH1, PH2, PH3**. They are the same tests: **M1 = PH1, M2 = PH2, M3 = PH3**. The thesis uses "PH". The CSV files, Python scripts and result plots still use "M". Both names appear below. The reference captures keep the names **R1 – R4**.

---

## 1. What is measured, and the main results

**Phase deviation** is how far the output phase moves when the message is applied. For analog phase modulation (PM) the output is

`s(t) = A cos( 2π fc t + kp · sin(2π fm t) )`

and `kp` (radians) is the phase deviation that the firmware is told to produce. **Phase accuracy** is how close the *measured* deviation is to the *commanded* `kp`.

Three independent methods are used, so that the result does not rest on one technique:

| Test | Name | What it checks | How |
|---|---|---|---|
| **PH1** (M1) | Commanded-phase sweep | Measured phase deviation vs. commanded, over kp = 0.25 … 4.0 rad | Hilbert transform gives the phase of each capture; a line is fitted to measured vs. commanded |
| **PH2** (M2) | Bessel null (J0) | The exact `kp` where the carrier disappears | At kp = 2.40483 rad the carrier amplitude is zero (J0 = 0). The scope amplitude is measured while `kp` is stepped around that point |
| **PH3** (M3) | BPSK phase reversal | A 180° phase step between two symbols | Phase of the "0" state and the "1" state, plus the transition time between them |

Four **reference captures** (R1 – R4) measure the things that would otherwise be mixed into the result: the true carrier frequency (R1), the delay and phase offset of the output filter (R2), the measurement noise floor (R3), and the envelope ripple (R4).

### Headline results (all reproducible from this folder)

| Quantity | Result |
|---|---|
| Carrier frequency, measured (R1) | **5000.4221 Hz** (+84.4 ppm). This value is used by all later scripts |
| Filter offset and delay (R2) | **−19.47°**, equivalent to **10.82 µs** delay |
| Noise floor (R3), one capture, no averaging | RMS 20.05 mV = 24.81 DAC codes → phase floor **1.18°**. With scope averaging: **0.30°** (×16), **0.15°** (×64) |
| Envelope ripple at kp = 1 (R4) | **2.82 %** (limit 5 %) |
| PH1 phase-sensitivity slope (30 captures) | **0.9946** (ideal 1.0000). Gain error **−0.54 %**, R² = 0.999988 |
| PH1 linearity | Largest deviation from the line **0.24 % of full scale** (full scale = 4 rad) |
| PH1 gain error per point | −0.03 % at kp ≤ 2, growing to **−0.64 %** at kp = 4 |
| PH2 first Bessel null | Fitted at **2.4073 rad** (theory 2.40483 rad), error **+0.10 %**, depth −39.8 dB at the fitted minimum |
| PH3 phase reversal | Separation **180.00° ± 0.16°** (3 repeats; ±180° is the same point on the circle) |
| PH3 transition time (10 % – 90 %) | **45.9 µs**, overshoot **1.9°** |
| PH3 independent check (zero-crossing timing) | **−179.1°** (Δt = −99.49 µs on a 200 µs period) |

---

## 2. Equipment and wiring

| Item | Use |
|---|---|
| STM32 Nucleo-F303RE | Waveform generator: DAC1_OUT1 on **PA4**, TIM6 trigger → circular DMA |
| Sallen-Key Butterworth low-pass filter (MCP6002, R = 5.6 kΩ, C = 2 nF / 1 nF, f_c ≈ 20.1 kHz, Q = 0.707) | Reconstruction filter. **All reported measurements are taken at its output** |
| Tektronix TBS1102C oscilloscope | 2 channels, 8-bit, **fixed 2000-point record** |
| 10× probes | CH1 on the filter output, CH2 on **PB0** (sync pulse). For R2, CH2 goes on raw **PA4** |
| USB flash drive | CSV export, one file per channel |
| STM32CubeIDE | Build and flash |
| Python 3 | numpy, scipy, matplotlib |

**Signal path:** `PA4 (DAC) → Sallen-Key LPF → scope CH1 → CSV → Python`.
**Sync path:** `PB0 → scope CH2` (a 40 µs pulse at every DMA buffer wrap, used as the trigger and, for PH3, as the time reference).

Use one common ground, keep the probe ground clip short, and power the laptop from battery or its charger alone (not USB + mains together), which otherwise adds visible oscillation to the noise measurements.

---

## 3. Firmware (STM32 code)

Both firmware files are in [`firmware/`](firmware). They differ from each other only in what they do with the DAC.

| File | Used for |
|---|---|
| [`firmware/phase_accuracy.c`](firmware/phase_accuracy.c) | **R1, R2, R4, PH1, PH2, PH3.** One file. You change a few `#define` lines (and `kp_idx`) for each test, rebuild and flash |
| [`firmware/R3_reference.c`](firmware/R3_reference.c) | **R3 only.** Holds the DAC at a fixed code (2048). No timer, no DMA, no sync pulse, so only noise is left |

Both files are the `main.c` of a CubeMX project. The peripheral set-up (DAC, DMA, GPIO, TIM6) comes from the CubeMX-generated `dac.c`, `dma.c`, `gpio.c`, `tim.c`. The header comment of each file lists what it defines and what it expects from the CubeMX files.

### 3.1 What `phase_accuracy.c` does

- Clock: HSE 8 MHz (bypass, from the ST-LINK) → PLL ×9 → **72 MHz**.
- Fills a buffer with `MID_CODE + AMP_CODES · cos(2π fc t + kp · sin(2π fm t))` and streams it to the DAC using **TIM6 trigger → circular DMA**. The CPU does nothing during output.
- DAC update rate is exactly **100 kSa/s** (TIM6: PSC = 0, ARR = 719).
- At every buffer wrap it raises **PB0** for `SYNC_PULSE_US` (40 µs), timed with the DWT cycle counter. This is the scope trigger.
- The VCP/UART output on this board is garbled, so the `kp` value is **not** read over UART. Choose it by setting `kp_idx` in the source and reflashing (recommended). Button B1 steps through the list, and LD2 blinks the index (one blink = index 0).

### 3.2 CONFIG block: what to set for each test

```c
#define MODE_ANALOG_PM   1          /* 1 = analog PM (PH1, PH2); 0 = BPSK (PH3) */
#define FC_HZ            5000.0f    /* carrier                                  */
#define FM_HZ            500.0f     /* message; 1000.0f for PH2 only            */
#define NSAMP            1000u      /* buffer length; 4000u for PH3             */
#define RS_SYM           250.0f     /* BPSK symbol rate                         */
#define AMP_CODES        1200       /* peak amplitude in DAC codes              */
#define MID_CODE         2048
#define SYNC_PULSE_US    40u

static const float kp_list[] = { 0.00f, 0.25f, 0.50f, 0.75f, 1.00f,
                                 1.50f, 2.00f, 2.50f, 3.00f, 3.50f, 4.00f };
static volatile uint8_t kp_idx = 0;   /* index into kp_list */
```

| Test | Firmware | Settings to change | Captures |
|---|---|---|---|
| **R1** carrier | `phase_accuracy.c` | `kp_idx = 0` (kp = 0) | 1 |
| **R2** filter offset | `phase_accuracy.c` | `kp_idx = 0` (kp = 0). CH1 = filter output, CH2 = raw PA4 | 2 files (`R2_CH1`, `R2_CH2`) |
| **R3** noise floor | `R3_reference.c` | none. Flash this one first, before the other firmware | 1 |
| **R4** envelope | `phase_accuracy.c` | `kp_idx = 4` (kp = 1.00) | 1 |
| **PH1** (M1) sweep | `phase_accuracy.c` | `kp_idx = 1 … 10` (kp = 0.25 … 4.00), one flash per point | 10 values × 3 repeats = **30** |
| **PH2** (M2) coarse | `phase_accuracy.c` | `FM_HZ 1000.0f`. Replace `kp_list` with `{2.20, 2.25, 2.30, 2.35, 2.40, 2.45, 2.50, 2.55, 2.60}` | **9** |
| **PH2** (M2) fine | `phase_accuracy.c` | `FM_HZ 1000.0f`. Replace `kp_list` with `{2.30, 2.32, … , 2.50}` (step 0.02) | **11** |
| **PH3** (M3) BPSK | `phase_accuracy.c` | `MODE_ANALOG_PM 0`, `NSAMP 4000u`, `RS_SYM 250.0f`. The DMA transfer count follows `NSAMP`. Bit pattern `{0,1,0,1,1,0,1,0,0,1}`, 4 ms per bit, so one buffer holds the 10 bits | 3 repeats + 1 zoomed |

> Only `phase_accuracy.c` for kp_idx = 0 is stored as-is in the repository. The PH2 and PH3 builds are this same file with the changes in the table above. That is what was flashed.

---

## 4. Oscilloscope settings

All captures use a **fixed 2000-point record** at 1 ms/div, which gives 8 µs per sample (125 kSa/s) and a 16 ms record (80 carrier cycles, 25 samples per cycle). The record length and sample interval are written in the CSV header, and the scripts check them.

| Setting | Value |
|---|---|
| Timebase | **1 ms/div** (PH3 transition zoom: **200 µs/div**, 1.6 µs per sample) |
| Vertical | **1 V/div**, left unchanged for every capture (do not rescale between R3 and the measurement captures, or the noise floor no longer applies) |
| Coupling | CH1 **AC** |
| Trigger | **CH2 (PB0 sync), falling edge, about 125 mV, Normal mode** for all modulated captures. A CH1 self-trigger is only valid for kp = 0 |
| Acquisition | Average ×16 for the modulated captures (×64 for PH2). R3 is a single **Sample** capture with Auto trigger (nothing to average to) |
| Export | USB flash drive, one CSV per channel. Hi-Res acquisition mode cannot be used (it is incompatible with the CSV transfer) |

Settings as recorded in the CSV headers:

| Files | Sample interval | Horizontal scale | Vertical scale | Probe setting in header |
|---|---|---|---|---|
| R1, R2_CH1, R2_CH2 | 8 µs | 1 ms/div | 1 V/div | 1× |
| R3, R4 | 8 µs | 1 ms/div | 1 V/div | 10× |
| PH1, PH2, PH3 | 8 µs | 1 ms/div | 1 V/div | 10× |
| PH3 zoomed transition | 1.6 µs | 0.2 ms/div | 1 V/div | 10× |

Before analysing any new CSV, check it with `diagnose.py` (Section 6.2).

---

## 5. Folder layout and file names

```
phase_accuracy/
├── csv/
│   ├── Reference/   R1.CSV  R2_CH1.CSV  R2_CH2.CSV  R3.CSV  R4.CSV
│   ├── M1/          m1_kp<K>_r<N>.CSV        (PH1: 10 kp values × 3 repeats)
│   ├── M2/
│   │   ├── J0 coarse pass/   m2_j0_coarse_<kp>.CSV   (9 files)
│   │   └── J0 fine pass/     m2_j0_fine_<kp>.CSV     (11 files)
│   └── M3/          m3_ch1_r<N>.CSV  m3_ch2_r<N>.CSV  (PH3: 3 repeats, CH1 + CH2)
│                    m3_transition.CSV  m3_transition_ch2.CSV  (zoomed)
├── firmware/        phase_accuracy.c  R3_reference.c
├── python/          all analysis scripts + stage.py
├── oscilloscope_images/   R1.png  m1_kp4.PNG  m2_kp2.6.PNG  m3_zoomed.PNG
└── result_plots/    15 PNG figures (Section 7)
```

**File-name rules**

- `kp` in PH1 names uses `p` for the decimal point: `m1_kp0p25_r1.CSV` = kp 0.25, repeat 1. `_r1 … _r3` are repeated captures.
- `m3_ch1_*` is the filter output, `m3_ch2_*` is the PB0 sync pulse, both from the same capture.
- The CSVs in this repository are sorted into sub-folders and have **upper-case `.CSV`** extensions. The scripts were written for one flat folder with lower-case `.csv` names (`m1_kp1p00_r1.csv`, `ref_carrier_kp0.csv`). `stage.py` converts between the two (next section).

---

## 6. Reproducing the results

### 6.0 Setup

```
pip install numpy scipy matplotlib
```

Tested with numpy 2.5, scipy 1.18, matplotlib 3.11 (older versions should also work).

### 6.1 Stage the files

Open a terminal in the `phase_accuracy` folder and run:

```
python python/stage.py
cd work
```

`stage.py` copies all scripts and CSVs into a new flat folder `work/`, renames files to the names the scripts expect (`R1.CSV` → `ref_carrier_kp0.csv`, and so on) and writes `m1_manifest.txt` and `m3_manifest.txt`. Nothing in the repository is changed. **All commands below are run inside `work/`.**

| Repository file | Name in `work/` |
|---|---|
| `csv/Reference/R1.CSV` | `ref_carrier_kp0.csv` |
| `csv/Reference/R2_CH1.CSV`, `R2_CH2.CSV`, `R3.CSV`, `R4.CSV` | same names |
| `csv/M1/m1_kp1p00_r1.CSV` | `m1_kp1p00_r1.csv` |
| `csv/M2/J0 coarse pass/m2_j0_coarse_2.40.CSV` | `m2_j0_coarse_2.40.csv` |
| `csv/M2/J0 fine pass/m2_j0_fine_2.40.CSV` | `m2_j0_fine_2.40.csv` |
| `csv/M3/m3_ch1_r1.CSV`, `m3_ch2_r1.CSV` | `m3_ch1_r1.csv`, `m3_ch2_r1.csv` |
| `csv/M3/m3_transition.CSV`, `m3_transition_ch2.CSV` | `m3_transition.csv`, `m3_transition_ch2.csv` |

### 6.2 Capture check (optional, for any new CSV)

```
python diagnose.py ref_carrier_kp0.csv
python diagnose.py m1_kp1p00_r1.csv --fc 5000 --interval 8e-6 --record 2000
```

Checks the record length, sample interval, clipping, and that the power is near 5 kHz. It ends with `PASS` or a warning list. Do not analyse a capture that does not PASS (R3 and the PH2 files near the null will show expected warnings; they are not sine waves at 5 kHz). `diagnose_csv_load.py <file>` shows how the CSV columns are read, for debugging a bad export.

### 6.3 Order of analysis

| Step | Script | Input CSV | Output |
|---|---|---|---|
| 1 | `measure_fc_ref.py` | R1 | carrier frequency, writes `fc_ref.json` |
| 2 | `check_r2.py` | R2 (CH1 + CH2) | filter offset and delay |
| 3 | `check_r3_noise.py` | R3 | noise floor and phase floor |
| 4 | `envelope_check.py` | R4 | envelope ripple |
| 5 | `check_m1_point.py` | one PH1 capture | β, gain error, RMS / peak phase error |
| 6 | `m1_sweep_summary.py` | all 30 PH1 captures (manifest) | slope, gain error, linearity |
| 7 | `quick_coarse_check.py` | PH2 coarse (9) + R1 | coarse null position |
| 8 | `m2_j0_fine_fit.py` | PH2 fine (11) + R1 | fitted null, error, depth |
| 9 | `check_m3_point.py` | one PH3 repeat (CH1 + CH2) | state phases, separation, ripple |
| 10 | `m3_sweep_summary.py` | 3 PH3 repeats (manifest) | mean ± spread over repeats |
| 11 | `m3_transition_detail.py` | zoomed PH3 | transition time, overshoot, figure |
| 12 | `m3_simple_demo.py` | PH3 CH1 | zero-crossing check, figure |
| 13 | `make_figures.py` | several | the 13 remaining figures |

All scripts after step 1 take the carrier frequency **5000.422114 Hz** from R1 as `--fc`. Use the value your own run prints.

### 6.4 Commands and expected output

**Step 1 — carrier frequency (R1)**

```
python measure_fc_ref.py ref_carrier_kp0.csv
```
```
zero crossing    : 5000.422114 Hz (+84.423 ppm vs nominal)    (77 crossings used)
phase slope      : 5000.378106 Hz
FFT parabolic    : 5000.404 Hz
implied DAC fs   : 100008.442 Sa/s (+84.42 ppm vs 100 kSa/s)
f_c_hat = 5000.422114 Hz  ->  written to fc_ref.json
```

The zero-crossing estimate is the one used. The script prints a warning that estimators 1 and 2 differ by more than 5 ppm (8.8 ppm here). The three estimates agree to within 9 ppm. A 9 ppm carrier error is a phase ramp of about 0.26° across the whole 16 ms record, which is below the 0.30° noise floor, so it does not change any phase result. The +84 ppm agrees with the HSE result of the frequency-accuracy metric.

**Step 2 — filter offset and delay (R2)**

```
python check_r2.py R2_CH1.CSV R2_CH2.CSV --fc 5000.422114
```
```
CH1 (LPF output) amplitude : 0.9698
CH2 (raw PA4)    amplitude : 0.9537
Constant phase offset       : -19.47 deg
Equivalent time delay       : -10.82 us
```

The prediction from the nominal filter values is about −20.7°; the measured −19.47° is 1.2° away, which is inside component tolerance. The DAC zero-order-hold delay is in both channels, so it cancels in this CH1-minus-CH2 difference: R2 measures the **filter only**.

**Step 3 — noise floor (R3)**

```
python check_r3_noise.py R3.CSV --codes-per-volt 1237.4
```
```
noise RMS : 20.049 mV  (24.809 codes)
sigma_phi = 1.1846 deg      Average x16: 0.2961 deg      Average x64: 0.1481 deg
```

`1237.4` = 1200 codes / 0.9698 V (from R2). Phase floor = noise (in codes) / amplitude (1200 codes). This is the **instrument-limited floor**: any RMS phase error at or below it cannot be told apart from noise.

**Step 4 — envelope ripple (R4)**

```
python envelope_check.py R4.CSV --fc 5000.422114
```
```
Ripple : 2.82 %  (limit: 5 %)      PASS
```

A flat envelope is a condition for the phase extraction to be valid.

**Step 5 — one PH1 point**

```
python check_m1_point.py m1_kp1p00_r1.csv --kp 1.0 --fc 5000.422114
```
```
Measured beta : 0.9997 rad       Gain error : -0.03 %
RMS phase error : 0.113 deg (noise floor 0.30 deg)       Peak : 0.310 deg
```

The script also prints a fitted delay and offset. These are **nuisance parameters**: they only line up the measured and ideal curves before the error is computed, and they are not compared with R2. An RMS error at or below the noise floor is reported as instrument-limited (quote 0.30°).

**Step 6 — full PH1 sweep**

```
python m1_sweep_summary.py m1_manifest.txt --fc 5000.422114
```
```
Slope (realised phase sensitivity) : 0.9946
Gain error                         : -0.54 %
Intercept                          : +0.0052 rad
R^2                                : 0.999988
Max linearity deviation            : 0.24 % of full scale
```

It first prints β and the gain error for all 30 captures (for example kp = 4.00 → β = 3.9746, −0.63 %). `m1_manifest.txt` lists "kp  filename" per line.

**Steps 7 – 8 — PH2 Bessel null**

```
python quick_coarse_check.py
python m2_j0_fine_fit.py --ref ref_carrier_kp0.csv
```
```
Coarse minimum: kp = 2.40  (-43.9 dB)

Fitted null location : kp = 2.4073 rad
Theoretical (J0)     : kp = 2.40483 rad
Gain error           : +0.10 %
Suppression at null  : -39.8 dB
```

The coarse script takes no arguments (it reads the `m2_j0_coarse_*.csv` files in the folder). The fine script has a `KP_LIST` at its top, which must match the list flashed to the board. `check_m2_point.py <file> --kp <kp> --fc 5000.422114 --ref ref_carrier_kp0.csv` checks a single PH2 capture.

**Steps 9 – 10 — PH3 BPSK**

```
python check_m3_point.py m3_ch1_r1.csv m3_ch2_r1.csv --fc 5000.422114
python m3_sweep_summary.py m3_manifest.txt --fc 5000.422114
```
```
Separation              : +180.00 +/- 0.16 deg  (ideal 180)
Mean phase, state 0     : +19.70 +/- 8.79 deg
Steady-state RMS error  : 19.72 +/- 8.77 deg
Envelope ripple (ss)    : 7.00 +/- 1.55 %
Guard interval +/-400 us, post-detection bandwidth 2500 Hz
```

CH2 (the PB0 pulse) gives the time origin for the ideal phase pattern. The **separation** between the two symbol phases is the accuracy result; the steady-state RMS is dominated by one constant phase offset that differs from capture to capture (the per-capture RMS values are 32.1°, 14.7° and 12.5°) and is not a measure of noise. Separation is averaged with circular statistics, because +180° and −180° are the same phase. `m3_manifest.txt` lists "ch1_file  ch2_file" per repeat.

**Step 11 — PH3 transition (zoomed capture)**

```
python m3_transition_detail.py m3_transition.csv --fc 5000.422114
```
```
Phase step               : +182.54 deg
Transition time (10-90%) : 45.92 us
Overshoot                : 1.90 deg
```

Writes `m3_transition_detail.png`. Look at the figure before trusting the numbers.

**Step 12 — zero-crossing check of the phase reversal**

```
python m3_simple_demo.py m3_ch1_r1.csv --fc 5000.422114
```
```
Delta t   = -99.49 us
Delta phi = -179.1 deg   (ideal: 180.0 deg)
```

This does not use a Hilbert transform. It compares the time of the carrier's zero-crossings before and after the symbol change: 99.49 µs on a 200 µs carrier period is 179.1°. It agrees with the Hilbert result to within 1°. Writes `m3_phase_reversal_demo.png`.

**Step 13 — figures**

```
python make_figures.py
```

Writes 13 PNG files to `work/figures/`. The file names at the top of the script (CONFIG block) already match the staged names. It needs `ref_carrier_kp0.csv`, `R3.CSV`, `m1_kp1p00_r1.csv`, `m1_manifest.txt`, the 11 fine-pass CSVs, and `m3_ch1_r1.csv` + `m3_ch2_r1.csv`.

---

## 7. Result plots: what each one shows

All plots are in [`result_plots/`](result_plots). "Fig." is the figure number in the thesis (Section 7.7). The two plots made by their own scripts are marked.

### PH1 (M1): commanded-phase sweep

| File | Fig. | Script | What it shows and what to look for |
|---|---|---|---|
| `pm_waveform.png` | 7.30 | `make_figures.py` | The 5 kHz output at kp = 1 rad (top) and the PB0 sync pulse (bottom). Spacing of the cycles is stretched and squeezed by the message, which is the phase modulation |
| `envelope.png` | 7.31 | `make_figures.py` | Envelope of the signal at kp = 1. It should be flat (ripple 0.54 % here), so the phase can be read without amplitude effects |
| `trajectory.png` | 7.32 | `make_figures.py` | Measured phase against time, overlaid on the ideal `kp·sin(2π fm t)` after removing delay and offset. The curves should lie on top of each other |
| `phase_error.png` | 7.33 | `make_figures.py` | Measured minus ideal phase against time (degrees). Should be small and without a pattern; a sine pattern at 500 Hz would mean a gain error |
| `histogram.png` | 7.34 | `make_figures.py` | Spread of those errors. A narrow, centred bell means random noise only |
| `m1_sweep.png` | 7.35 | `make_figures.py` | **Key plot.** Measured β against commanded kp (30 points), with the fitted line and the ideal slope-1 line. Slope 0.9946, R² = 0.999988 |
| `m1_residuals.png` | 7.36 | `make_figures.py` | Distance of each point from the fitted line (mrad). A smooth curve means a small systematic effect; random scatter means noise |
| `gain_error.png` | 7.37 | `make_figures.py` | Gain error (%) against kp. About −0.03 % up to kp = 2, then growing to −0.64 % at kp = 4. The error grows with the size of the phase swing |

### PH2 (M2): Bessel null

| File | Fig. | Script | What it shows and what to look for |
|---|---|---|---|
| `m2_fine.png` | 7.38 | `make_figures.py` | **Key plot.** Carrier amplitude (dB, relative to unmodulated) against kp from 2.30 to 2.50, with a parabola fitted to the dip. The fitted minimum is at 2.4073 rad, compared with the theoretical 2.40483 rad |
| `spectra.png` | 7.39 | `make_figures.py` | Spectrum at kp = 0 and at the null. At the null the 5 kHz carrier is gone while the sidebands stay |

### PH3 (M3): BPSK phase reversal

| File | Fig. | Script | What it shows and what to look for |
|---|---|---|---|
| `bpsk_waveform.png` | 7.40 | `make_figures.py` | BPSK output (top) and PB0 sync (bottom). The carrier flips phase at each bit change |
| `bpsk_trajectory.png` | 7.41 | `make_figures.py` | Measured phase against the ideal two-level pattern. Grey bands are the ±400 µs guard intervals around each transition, which are left out of the steady-state statistics |
| `m3_transition_detail.png` | 7.42 | `m3_transition_detail.py` | Zoomed look at one 180° step. Marks the 10 % and 90 % points (45.9 µs) and the small overshoot (1.9°) |
| `m3_phase_reversal_demo.png` | 7.43 | `m3_simple_demo.py` | Zero-crossings of the carrier before and after the symbol change. The time shift of 99.49 µs corresponds to 179.1° |

### Noise floor (R3)

| File | Fig. | Script | What it shows and what to look for |
|---|---|---|---|
| `noise.png` | 7.44 | `make_figures.py` | Left: the R3 capture with the DAC held at mid-scale (RMS 20.05 mV). Right: its spectrum. It should look like flat noise with no strong single tone |

### Oscilloscope images

[`oscilloscope_images/`](oscilloscope_images) holds screen captures from the instrument:

| File | Thesis fig. | Content |
|---|---|---|
| `R1.png` | (not numbered) | Unmodulated 5 kHz reference capture (R1) |
| `m1_kp4.PNG` | 7.37A | PH1 at kp = 4.0 rad |
| `m2_kp2.6.PNG` | 7.39A | PH2 coarse pass at kp = 2.6 rad |
| `m3_zoomed.PNG` | 7.41A | PH3 zoomed phase transition |

---

## 8. Method details

- **Phase extraction (PH1, PH3).** The mean is removed and the Hilbert transform gives the complex (analytic) signal. It is multiplied by `exp(−j2π fc t)` using the R1 carrier, filtered with a zero-phase 4th-order Butterworth low-pass (post-detection bandwidth **3 kHz** for analog PM, **2.5 kHz** for BPSK, **15 kHz** for the zoomed transition) and the angle is taken with unwrapping. The first and last **10 carrier cycles** are discarded to remove filter edge effects.
- **Delay and offset (PH1).** Trigger timing and filter delay are not known exactly, so they are fitted as nuisance parameters before the error against the ideal `kp·sin(2π fm t)` is computed. β is the fitted amplitude of the measured phase at fm. The slope across the whole sweep is the "realised phase sensitivity".
- **Bessel null (PH2).** The carrier amplitude of a PM signal is proportional to J0(kp). The amplitude at 5 kHz is read from a Blackman-windowed FFT (fm = 1 kHz, so the sidebands are 1 kHz apart and do not leak into the carrier bin), compared with the R1 amplitude, and a parabola is fitted to the dB values around the minimum. The minimum position is much less sensitive to amplitude calibration than a direct amplitude measurement.
- **BPSK (PH3).** The PB0 sync pulse gives the time of the buffer start, and from it the ideal phase pattern for the bit sequence `0 1 0 1 1 0 1 0 0 1` at 250 symbols/s. The mean phase of each state is found with the guard interval removed, and the separation between the two states is computed with circular statistics. The independent zero-crossing check does not use any filtering.
- **Instrument limit.** R3 shows the phase floor is 1.18° (single capture), 0.30° (×16) and 0.15° (×64). RMS phase errors at or below this are instrument-limited.

---

## 9. Known limitations and data notes

1. **PH1 RMS and peak errors per point.** The sweep summary reports β and gain error for all 30 captures, but not RMS and peak phase error. The thesis quotes these for kp = 0.25 and kp = 1.0. `check_m1_point.py` gives them for any single capture.
2. **J1 (second Bessel null)** was not measured. Only the first null (J0) was used.
3. **PH3 steady-state RMS** (about 20°) is a constant per-capture phase offset, not noise. Use the separation (180.00° ± 0.16°) as the accuracy result. Two of the three repeats have steady-state envelope ripple above 5 % (8.67 % and 7.38 %); the third is 4.93 %.
4. **PH3 repeat spread** is the standard deviation of three captures, which is small for a confidence statement.
5. **R1 and R2 headers** record probe setting 1×, while all other files record 10×. The recorded amplitudes (0.976 V and 0.970 V) are correct for the 1 V/div setting, so the values are in real volts. Phase results do not depend on this.
6. **`check_r2.py`** still prints a generic expectation of about −30° (filter plus zero-order hold). The measured value of −19.47° is the filter alone (the zero-order-hold delay cancels between CH1 and CH2). The "PASS" message is correct; the "near −30°" text is outdated.
7. **Frequency reference.** The scope timebase is ±20 ppm. A ±20 ppm frequency error is a phase ramp of about 0.06° across one 16 ms record, far below the noise floor, so it does not limit any result here.
8. **PH2** depth (−39.8 dB) is the fitted parabola minimum. The deepest single capture is −43.8 dB (kp = 2.40), and the depth is limited by the 8-bit scope and by the step size of 0.02 rad.
9. **`check_m1_point.py`** prints a comment that PH1 captures are "self-triggered". The captures were triggered on the PB0 sync pulse, but this does not affect the result, because delay and offset are fitted for every capture.

---

## 10. Quick start (no hardware)

```
git clone https://github.com/moyeenurrahman/STM32_phase_modulation_research
cd STM32_phase_modulation_research/phase_accuracy
pip install numpy scipy matplotlib
python python/stage.py
cd work

python measure_fc_ref.py ref_carrier_kp0.csv                       # R1: carrier
python check_r2.py R2_CH1.CSV R2_CH2.CSV --fc 5000.422114          # R2: filter
python check_r3_noise.py R3.CSV --codes-per-volt 1237.4            # R3: noise floor
python envelope_check.py R4.CSV --fc 5000.422114                   # R4: ripple
python m1_sweep_summary.py m1_manifest.txt --fc 5000.422114        # PH1 result
python quick_coarse_check.py                                       # PH2 coarse
python m2_j0_fine_fit.py --ref ref_carrier_kp0.csv                 # PH2 result
python m3_sweep_summary.py m3_manifest.txt --fc 5000.422114        # PH3 result
python m3_transition_detail.py m3_transition.csv --fc 5000.422114  # PH3 transition
python m3_simple_demo.py m3_ch1_r1.csv --fc 5000.422114            # PH3 check
python make_figures.py                                             # remaining figures
```

Every number in Section 1 should appear in the output.
