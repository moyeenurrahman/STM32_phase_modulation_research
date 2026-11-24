# SDR Extension: Analysis and Reproduction Guide

Part of the thesis *"Implementation and Performance Analysis of a Baseband Phase Modulation (PM) System with STM32 for Low-Power Wireless Communication"* (BSc EEE, RUET). This folder holds the work of **Chapter 8 (Software-Defined Radio extension)**.

This guide explains three things:

1. **QPSK transmission.** How the STM32 generated a QPSK signal (NRZ and raised-cosine pulse shapes) and how it was captured.
2. **Software-defined receiver.** How a Python receiver demodulates the captured QPSK and BPSK signals and measures EVM. This includes an EVM-vs-SNR test made by adding noise to a real capture.
3. **BER simulation.** How a Monte Carlo simulation gives the BPSK bit error rate against Eb/N0 and checks it against theory.

For each result, it tells you which firmware was flashed, how the oscilloscope was set, which CSV files were recorded, which Python command was run, and which plot shows what. **Every number and plot can be reproduced from the files in this folder, without any hardware.** All commands in Section 6 were run on a clean copy of this folder and gave the same numbers as the stored result files.

> **Why "SDR"?** The transmitter is already software-defined: the waveform is a table of numbers computed in firmware and played through a plain 12-bit DAC. Changing from PM to BPSK to QPSK is a firmware change, not a hardware change. The receiver here is also software-defined: it takes only the sampled waveform (a scope CSV, standing in for an ADC) and does all demodulation in Python. Everything stays at baseband; there is no RF stage.

---

## 1. What is measured, and the main results

| Question | Data | How it is answered |
|---|---|---|
| Can the same STM32 DAC/timer/DMA system send **QPSK**? How clean is it? | QPSK NRZ, 3 captures | Software receiver → constellation, **EVM**, MER, phase error, bit errors |
| Does **raised-cosine** pulse shaping improve the signal? | QPSK RC, 3 captures | Same receiver; EVM, spectrum and eye diagram compared with NRZ |
| Does the same receiver work on a **different signal it was not written for**? | BPSK, 5 captures (from the spectral-analysis test) | Same receiver chain with BPSK settings → EVM, bit errors, carrier frequency |
| How does EVM change with **noise**? | One clean BPSK capture + synthetic noise | EVM-vs-SNR sweep, compared with theory |
| What is the true **BPSK bit error rate** against Eb/N0? | No CSV (simulation) | Monte Carlo simulation, checked against `Q(√(2·Eb/N0))` |

### Headline results (all reproducible from this folder)

| Quantity | Result | Thesis ref. |
|---|---|---|
| QPSK NRZ EVM (RMS), 3 captures, 16 symbols | **1.636 %** (MER 35.72 dB, phase error 0.809° RMS), 0 of 32 bit errors | Table 8.3 |
| QPSK raised-cosine EVM (RMS), 3 captures, 16 symbols | **0.514 %** (MER 45.78 dB, phase error 0.195° RMS), 0 of 32 bit errors. Below the ≈0.6 % scope floor, so reported as an **upper bound** | Table 8.3 |
| QPSK carrier found by the receiver | 1000.4441 Hz (+444 ppm from the nominal 1 kHz) | Section 8.3 |
| RC vs NRZ spectrum | RC sidelobes fall to about −60 dB just outside fc ± Rs; NRZ sidelobes stay at about −20 to −40 dB | Fig. 8.5 |
| BPSK EVM (RMS), 5 captures, 75 symbols | **1.972 %** (MER 34.10 dB, phase error 0.321° RMS), 0 of 75 bit errors | Table 8.4 |
| BPSK carrier found by the receiver | 5014.9591 Hz = **+2992 ppm**. The HSI clock error measured in the frequency-accuracy test was +3126 ppm: same sign and size | Section 8.4 |
| EVM-vs-SNR sweep | Within about 2.5–3 dB of matched-filter theory at low SNR; flattens near 2.4 % at 20 dB, close to the 1.972 % clean-capture floor | Fig. 8.10 |
| BER simulation check at Eb/N0 = 0 dB | 0.0790 simulated vs 0.0786 theory (0.4 %) → **PASS** | Section 8.5 |
| BER simulation, 0 – 10.5 dB | 22 points, each with ≥ 100 errors; simulated/theory ratio 0.935 – 1.159 across 7 decades of BER | Table 8.5, Fig. 8.11 |

---

## 2. Equipment and wiring

| Item | Use |
|---|---|
| STM32 Nucleo-F303RE | Waveform generator: DAC1_OUT1 on **PA4**, TIM6 TRGO → circular DMA |
| Sallen-Key Butterworth low-pass filter (MCP6002 at 4.05 V, R = 5.6 kΩ, C1 = 2 nF, C2 = 1 nF, f_c ≈ 20.1 kHz) | Reconstruction filter, the same one used in every other test |
| Tektronix TBS1102C oscilloscope | CH1 only, 8-bit, fixed 2000-point record |
| STM32CubeIDE | Build and flash |
| Python 3 | numpy, scipy, matplotlib |

**Test point.** All SDR captures are at **Point B** (filter output). PA4 → filter input; filter output → scope CH1; common ground; short probe ground lead. Only PA4 is used. PA5 is not used (see Section 3.1).

---

## 3. Firmware (STM32 code)

| File | Signal | Used for |
|---|---|---|
| [`firmware/QPSK.c`](firmware/QPSK.c) | QPSK, passband, single DAC channel | `csv/QPSK_NRZ/` and `csv/QPSK_raised_cosine/` |
| [`firmware/BPSK_from_spectral_analysis.c`](firmware/BPSK_from_spectral_analysis.c) | BPSK, hard 0°/180° switching | `csv/BPSK_from_spectral_analysis/` |

Both files are the `main.c` of a CubeMX project. Both use **HSI 8 MHz with the PLL off**; the TIM6 reload value is computed at run time from the real bus clock, so the DAC rate is exactly 100 kSa/s.

### 3.1 `QPSK.c`

The firmware builds one buffer of QPSK samples at start-up and plays it forever through circular DMA. The CPU does nothing while the signal runs.

| Parameter | Value |
|---|---|
| DAC rate `fs` | 100 kSa/s |
| Carrier `fc` | 1 kHz |
| Symbol rate `Rs` | 200 Bd (2 bits/symbol → 400 bit/s) |
| Samples per symbol | 500 |
| Symbol pattern | `{0, 1, 3, 2, 0, 3, 1, 2}` (8 symbols, repeats) |
| Buffer | 4000 samples = 40 ms = exactly 40 carrier cycles (no jump at the DMA wrap) |
| Amplitude | Peak scaled to 2048 ± 800 codes in both pulse shapes |

The signal is computed as `s[n] = I·cos(2π·fc·n/fs) − Q·sin(2π·fc·n/fs)` with Gray-coded symbols:

| Symbol | Bits | I | Q | Phase |
|---|---|---|---|---|
| 0 | 00 | +1 | +1 | 45° |
| 1 | 01 | −1 | +1 | 135° |
| 2 | 10 | +1 | −1 | 315° |
| 3 | 11 | −1 | −1 | 225° |

**Lines to change** (then Clean → Build → Flash):

```c
#define QPSK_MODE     0    // 0 = passband on PA4. Always 0 for this data.
#define PULSE_SHAPE   0    // 0 = NRZ, 1 = raised cosine (beta = 0.35, +/-3 symbols)
```

| Build | `QPSK_MODE` | `PULSE_SHAPE` | CSV folder |
|---|---|---|---|
| QPSK NRZ | 0 | **0** | `csv/QPSK_NRZ/` |
| QPSK raised-cosine | 0 | **1** | `csv/QPSK_raised_cosine/` |

The file in this repository is saved with `PULSE_SHAPE 0`. Set it to 1 to remake the RC signal. With `PULSE_SHAPE 1` the buffer takes a few seconds to compute at start-up at 8 MHz; this is normal.

**Mode 1 (baseband I on PA4, Q on PA5) was not used.** On the Nucleo board PA5 also drives the green LED LD2, which clamps the output above about 1.8 V. The STM32F303RE has only one DAC (DAC1, outputs PA4 and PA5); there is no DAC2/PA6 on this chip. The comment near line 44 of `QPSK.c` that mentions PA6 (DAC2_OUT1) is out of date and does not apply to this board.

### 3.2 `BPSK_from_spectral_analysis.c`

This is the spectral-analysis firmware ([`../spectrum/firmware/spectral_analysis.c`](../spectrum/firmware/spectral_analysis.c)) with only one line changed:

```c
#define SPEC_MODE   MODE_BPSK_ABRUPT     // LED flashes 6 times at start-up
```

| Parameter | Value |
|---|---|
| `fs` / `fc` | 100 kSa/s / 5 kHz |
| Symbol rate | 500 bit/s (200 samples per symbol) |
| Buffer | 1000 samples = 10 ms = 5 symbols |
| Bit pattern | `BPSK_PATTERN 0x16` read from bit 0 upward → **0, 1, 1, 0, 1** (repeats every 10 ms) |
| Amplitude | 2048 ± 1200 codes |

These are the same captures as spectral case **S6N**. The 5 CSV files are byte-identical (MD5) to [`../spectrum/csv/S6N/`](../spectrum/csv/S6N). They were reused here; no new BPSK capture was made for the SDR work.

---

## 4. Oscilloscope settings

| Setting | Value |
|---|---|
| Channel | CH1 only, Point B, DC coupling |
| Vertical | 500 mV/div |
| Timebase | **2 ms/div** |
| Sample rate / interval | 62.5 kSa/s / 16 µs |
| Record | 2000 points ≈ 32 ms |
| Trigger | CH1 rising edge, level near the middle of the swing (about 1.1–1.2 V) |
| Acquisition | **Single** for every capture |
| Captures | QPSK NRZ: 3, QPSK RC: 3, BPSK: 5 |

**About the file names.** The QPSK files are named `..._2p5msdiv_...` because 2.5 ms/div was planned. The TBS1102C only offers 2 ms and 4 ms in that range, so **2 ms/div was actually used**. The real sample rate comes from the CSV time column (16 µs), and the scripts use that, not the file name.

**What this means for QPSK.** One QPSK pattern is 40 ms, but a 2 ms/div record is about 32 ms. So each capture holds about 80 % of one pattern (5–6 usable symbols). This is why 3 captures are pooled (Section 8.2). The BPSK pattern is only 10 ms, so each BPSK capture holds about three repeats (15 usable symbols).

**Trigger behaviour.** With a phase-switched signal, a level trigger sometimes catches a 0° symbol and sometimes a 180° symbol, so the live trace jumps. This is normal. Each Single capture is still valid; the receiver finds the phase on its own.

**CSV format.** Each file has the Tektronix header (`Record Length`, `Sample Interval`, …) in the left columns and the time and voltage in the right columns. The scripts find the data block automatically. They use **1986 of the 2000 points** (the first rows share their line with header text and are skipped). This does not change the results.

---

## 5. Folder layout and which files go together

```
SDR/
├── csv/
│   ├── QPSK_NRZ/                       qpsk_mode0_nrz_2p5msdiv_cap1..3.csv
│   ├── QPSK_raised_cosine/             qpsk_mode0_rc_2p5msdiv_cap1..3.csv
│   └── BPSK_from_spectral_analysis/    1.csv .. 5.csv   (= spectrum/csv/S6N)
├── firmware/
│   ├── QPSK.c
│   └── BPSK_from_spectral_analysis.c
├── python/
│   ├── qpsk_analysis.py        QPSK software receiver
│   ├── bpsk_analysis.py        BPSK software receiver + EVM-vs-SNR sweep
│   ├── ber_simulation.py       BPSK BER Monte Carlo simulation (no CSV)
│   └── make_nrz_rc_figures.py  builds the three NRZ-vs-RC comparison figures
├── oscilloscope_images/        scope screenshots (Section 7.2)
├── result_plots/               plots used in the thesis (Section 7.1)
├── results_from_python/        text output of each run (Section 6)
├── reproduced_plots/           created when you run the commands below
└── SDR_analysis_procedure.md   this guide
```

| Firmware | CSV folder | Script | Run name (`--prefix`) | Result text file |
|---|---|---|---|---|
| `QPSK.c`, `PULSE_SHAPE 0` | `QPSK_NRZ/` | `qpsk_analysis.py` | `qpsk_nrz` | `qpsk_nrz_results.txt` |
| `QPSK.c`, `PULSE_SHAPE 1` | `QPSK_raised_cosine/` | `qpsk_analysis.py` | `qpsk_rc` | `qpsk_rc_results.txt` |
| `BPSK_from_spectral_analysis.c` | `BPSK_from_spectral_analysis/` | `bpsk_analysis.py` (carrier estimated) | `bpsk_abrupt` | `bpsk_abrupt_results.txt` |
| same | same | `bpsk_analysis.py` (carrier fixed + SNR sweep) | `bpsk_abrupt_snr` | `bpsk_abrupt_snr_results.txt`, `bpsk_abrupt_snr_snr_sweep.txt` |
| none | none | `ber_simulation.py` | `ber_bpsk` | `ber_bpsk_results.txt` |
| none (reads the QPSK outputs) | none | `make_nrz_rc_figures.py` | — | — |

---

## 6. Reproducing the results

### 6.0 Setup

```
git clone https://github.com/moyeenurrahman/STM32_phase_modulation_research
cd STM32_phase_modulation_research/SDR
pip install numpy scipy matplotlib
mkdir reproduced_plots
```

Tested with Python 3.13, numpy 2.5.3, scipy 1.18.1 and matplotlib 3.11.2. **Run every command from inside the `SDR` folder.** The scripts write their outputs to the path given by `--prefix`, so `--prefix reproduced_plots/...` keeps your outputs apart from the stored ones.

### 6.1 QPSK, NRZ

```
python python/qpsk_analysis.py --csv csv/QPSK_NRZ/qpsk_mode0_nrz_2p5msdiv_cap1.csv csv/QPSK_NRZ/qpsk_mode0_nrz_2p5msdiv_cap2.csv csv/QPSK_NRZ/qpsk_mode0_nrz_2p5msdiv_cap3.csv --mode passband --shape nrz --prefix reproduced_plots/qpsk_nrz
```

Expected printout:

```
qpsk_mode0_nrz_2p5msdiv_cap1.csv: rotation=270 deg, symbols=5, bit errors=0 of 10
qpsk_mode0_nrz_2p5msdiv_cap2.csv: rotation= 90 deg, symbols=5, bit errors=0 of 10
qpsk_mode0_nrz_2p5msdiv_cap3.csv: rotation=270 deg, symbols=6, bit errors=0 of 12
Carrier estimate      : 1000.4441 Hz (+444 ppm vs nominal)
Symbols analysed      : 16  (32 bits)
EVM (RMS)             : 1.636 %
EVM (peak)            : 3.049 %
MER                   : 35.72 dB
Phase error, RMS      : 0.809 deg
I/Q gain imbalance    : -0.051 dB
Bit errors            : 0 of 32
```

### 6.2 QPSK, raised-cosine

**First, see why the carrier must be fixed** (optional):

```
python python/qpsk_analysis.py --csv csv/QPSK_raised_cosine/qpsk_mode0_rc_2p5msdiv_cap1.csv csv/QPSK_raised_cosine/qpsk_mode0_rc_2p5msdiv_cap2.csv csv/QPSK_raised_cosine/qpsk_mode0_rc_2p5msdiv_cap3.csv --mode passband --shape rc --no-plots --prefix reproduced_plots/qpsk_rc_auto
```

This prints `WARNING: carrier search found no clearly dominant peak`, finds the carrier at 950.0980 Hz (wrong, at the edge of the search window), and gives about 96 % EVM with 7 bit errors. The raised-cosine signal has a changing envelope, so the 4th-power carrier line is weak in a 32 ms record. The script detects this and warns.

**Then run it properly**, using the carrier and symbol rate already found from the NRZ capture (same board, same clock, only the pulse shape differs):

```
python python/qpsk_analysis.py --csv csv/QPSK_raised_cosine/qpsk_mode0_rc_2p5msdiv_cap1.csv csv/QPSK_raised_cosine/qpsk_mode0_rc_2p5msdiv_cap2.csv csv/QPSK_raised_cosine/qpsk_mode0_rc_2p5msdiv_cap3.csv --mode passband --shape rc --fc 1000.4441 --fc-fixed --rs 200.0888 --prefix reproduced_plots/qpsk_rc
```

Use both `--fc-fixed` and `--rs 200.0888`. With `--fc-fixed` the script does not rescale the symbol rate, so the rate must be given already corrected (200 × 1000.4441/1000).

Expected printout:

```
qpsk_mode0_rc_2p5msdiv_cap1.csv: rotation=270 deg, symbols=5, bit errors=0 of 10
qpsk_mode0_rc_2p5msdiv_cap2.csv: rotation=270 deg, symbols=6, bit errors=0 of 12
qpsk_mode0_rc_2p5msdiv_cap3.csv: rotation=  0 deg, symbols=5, bit errors=0 of 10
EVM (RMS)             : 0.514 %
EVM (peak)            : 0.993 %
MER                   : 45.78 dB
Phase error, RMS      : 0.195 deg
I/Q gain imbalance    : +0.077 dB
Bit errors            : 0 of 32
NOTE: EVM is below the 8-bit instrument floor ... Report this as an upper bound
```

### 6.3 NRZ vs RC comparison figures

Run after 6.1 and 6.2:

```
python python/make_nrz_rc_figures.py
```

It reads `qpsk_nrz_results.txt`, `qpsk_rc_results.txt` and the four spectrum/eye PNGs from `reproduced_plots/` (or from the stored folders if you skipped 6.1–6.2) and writes `nrz_rc_comparison_bars.png`, `nrz_rc_spectrum_comparison.png` and `nrz_rc_eye_comparison.png`.

### 6.4 BPSK receiver (carrier estimated)

```
python python/bpsk_analysis.py --csv csv/BPSK_from_spectral_analysis/1.csv csv/BPSK_from_spectral_analysis/2.csv csv/BPSK_from_spectral_analysis/3.csv csv/BPSK_from_spectral_analysis/4.csv csv/BPSK_from_spectral_analysis/5.csv --fc 5000 --rs 500 --pattern 0,1,1,0,1 --shape nrz --prefix reproduced_plots/bpsk_abrupt
```

`--pattern 0,1,1,0,1` is required. It is the firmware's bit pattern (Section 3.2). Without it the script uses a placeholder and warns that the result cannot be trusted.

Expected printout:

```
1.csv: rotation=180 deg, symbols=15, bit errors=0 of 15
2.csv: rotation=  0 deg, symbols=15, bit errors=0 of 15
3.csv: rotation=180 deg, symbols=15, bit errors=0 of 15
4.csv: rotation=180 deg, symbols=15, bit errors=0 of 15
5.csv: rotation=  0 deg, symbols=15, bit errors=0 of 15
Carrier estimate      : 5014.9591 Hz (+2992 ppm vs nominal)
Symbols analysed      : 75  (75 bits)
EVM (RMS)             : 1.967 %
MER                   : 34.12 dB
Bit errors            : 0 of 75
```

The **+2992 ppm** carrier offset is the HSI clock error. Compare it with the frequency-accuracy result (+3126 ppm, [`../frequency_accuracy/`](../frequency_accuracy)).

### 6.5 BPSK receiver (carrier fixed) and EVM-vs-SNR sweep

This is the run reported in the thesis (Table 8.4, Fig. 8.10). The carrier and symbol rate from 6.4 are fixed, because at low SNR the automatic carrier search becomes unreliable.

```
python python/bpsk_analysis.py --csv csv/BPSK_from_spectral_analysis/1.csv csv/BPSK_from_spectral_analysis/2.csv csv/BPSK_from_spectral_analysis/3.csv csv/BPSK_from_spectral_analysis/4.csv csv/BPSK_from_spectral_analysis/5.csv --fc 5014.9591 --rs 501.4959 --fc-fixed --pattern 0,1,1,0,1 --shape nrz --prefix reproduced_plots/bpsk_abrupt_snr --snr-sweep --snr-min -5 --snr-max 20 --snr-step 1 --snr-realizations 50
```

Runtime: about 30 seconds. Expected printout (first part):

```
EVM (RMS)             : 1.972 %
EVM (peak)            : 3.328 %
MER                   : 34.10 dB
Phase error, RMS      : 0.321 deg
Residual quadrature   : 0.569 %
Bit errors            : 0 of 75
```

The 1.972 % here and 1.967 % in 6.4 differ slightly because in 6.4 each capture gets its own carrier estimate, while here one fixed carrier is used for all five captures. The thesis uses 1.972 % because it belongs to the same run as the sweep.

Then the sweep table (selected rows; the full table is in `bpsk_abrupt_snr_snr_sweep.txt`):

```
  SNR (dB)   EVM mean (%)   EVM std (%)
      -5.0         21.453         2.642
       0.0         12.428         1.532
      10.0          4.246         0.566
      20.0          2.415         0.194
```

The noise is seeded (`--snr-seed 0`), so the table is identical every time.

### 6.6 BPSK BER-vs-Eb/N0 simulation

```
python python/ber_simulation.py --ebno-min 0 --ebno-max 12 --ebno-step 0.5 --min-errors 100 --prefix reproduced_plots/ber_bpsk
```

Runtime: about 10 seconds. Expected printout:

```
Eb/N0 = 0 dB: simulated BER = 0.0790 (15791 errors / 200000 bits), theory = 0.0786
Relative difference: 0.4%  ->  PASS
 Eb/N0(dB)    BER (sim)   errors         bits  BER (theory)    ratio
       0.0    7.853e-02    15707       200000     7.865e-02    0.999
       ...
      10.5    1.256e-06      100     79600000     1.084e-06    1.159
      11.0    2.700e-07       27    100000000     2.613e-07    1.033  (bit ceiling hit)
      11.5    1.100e-07       11    100000000     5.329e-08    2.064  (bit ceiling hit)
      12.0    1.000e-08        1    100000000     9.006e-09    1.110  (bit ceiling hit)
```

The simulation is seeded (`--seed 0`), so the table is identical to `results_from_python/ber_bpsk_results.txt`. If the 0 dB check fails, the script stops and runs no sweep.

### 6.7 Output files

| Command | Writes to `reproduced_plots/` | Same as stored file |
|---|---|---|
| 6.1 | `qpsk_nrz_{waveform,spectrum,iq_tracks,constellation,eye}.png`, `qpsk_nrz_results.txt` | `result_plots/qpsk_nrz_*.png`, `results_from_python/qpsk_nrz_results.txt` |
| 6.2 | `qpsk_rc_*` (same five plots + results) | `result_plots/qpsk_rc_*.png`, `results_from_python/qpsk_rc_results.txt` |
| 6.3 | `nrz_rc_comparison_bars.png`, `nrz_rc_spectrum_comparison.png`, `nrz_rc_eye_comparison.png` | `result_plots/nrz_rc_*.png` |
| 6.4 | `bpsk_abrupt_*` (five plots + results) | `result_plots/bpsk_abrupt_*.png`, `results_from_python/bpsk_abrupt_results.txt` |
| 6.5 | `bpsk_abrupt_snr_*` (five plots, `_snr_sweep.png`, `_results.txt`, `_snr_sweep.txt`) | `result_plots/bpsk_abrupt_snr_*.png`, `results_from_python/bpsk_abrupt_snr_*.txt` |
| 6.6 | `ber_bpsk_curve.png`, `ber_bpsk_results.txt` | `result_plots/ber_bpsk_curve.png`, `results_from_python/ber_bpsk_results.txt` |

(`bpsk_abrupt_iq_tracks.png` from 6.4 is not stored in `result_plots/`; the thesis uses the 6.5 version.)

---

## 7. What each plot and screenshot shows

### 7.1 Result plots ([`result_plots/`](result_plots))

**QPSK** (`qpsk_nrz_*` from 6.1, `qpsk_rc_*` from 6.2)

| File | Thesis fig. | What it shows and what to look for |
|---|---|---|
| `qpsk_nrz_waveform.png` | 8.2(a) | Captured Point B voltage, about 1.2–2.45 V. Narrow notches where the phase jumps: deeper notches are 180° jumps, shallow ones are 90° jumps |
| `qpsk_rc_waveform.png` | 8.2(b) | Same for raised-cosine. The envelope shrinks smoothly toward the middle level (about 1.8 V) at 180° transitions instead of jumping |
| `qpsk_nrz_iq_tracks.png`, `qpsk_rc_iq_tracks.png` | 8.3 | The baseband I (blue) and Q (red) after the receiver mixes the signal down. Circles mark the chosen sampling instants; each circle sits on its curve |
| `qpsk_nrz_constellation.png` | 8.4(a) | Four tight clusters on the four Gray-coded points (crosses); EVM = 1.64 % in the title. Short red lines are the error vectors |
| `qpsk_rc_constellation.png` | 8.4(b) | Same, EVM = 0.51 % |
| `qpsk_nrz_spectrum.png`, `qpsk_rc_spectrum.png` | 8.5 (inputs) | Blackman-Harris spectrum, 0–5 kHz. Red dashed line: recovered carrier (1000.44 Hz). Grey dotted lines: fc ± Rs |
| `qpsk_nrz_eye.png`, `qpsk_rc_eye.png` | 8.6 (inputs) | In-phase eye diagram over two symbols. Each comes from one capture only, so it has few traces |
| `nrz_rc_spectrum_comparison.png` | **8.5** | **Key plot.** NRZ and RC spectra side by side. RC sidelobes fall to about −60 dB right outside fc ± Rs; NRZ sidelobes stay at about −20 to −40 dB |
| `nrz_rc_eye_comparison.png` | **8.6** | NRZ and RC eye diagrams side by side. RC traces are smoother through the crossing |
| `nrz_rc_comparison_bars.png` | **8.7** | Bars of EVM (RMS and peak), MER and RMS phase error. Red dashed line: ≈0.6 % scope EVM floor. NRZ is above it (measured value); RC RMS is below it (upper bound) |

**BPSK** (`bpsk_abrupt_snr_*` from 6.5; `bpsk_abrupt_*` from 6.4 are the same plots from the carrier-estimated run)

| File | Thesis fig. | What it shows and what to look for |
|---|---|---|
| `bpsk_abrupt_snr_constellation.png` | **8.8** | Two tight clusters at ±1 on the I axis; EVM = 1.97 % |
| `bpsk_abrupt_snr_waveform.png` | 8.9 | 5 kHz carrier with a notch at every phase reversal. No notch between two equal bits, as expected for the pattern 0,1,1,0,1 |
| `bpsk_abrupt_snr_spectrum.png` | 8.9 | sinc-shaped main lobe around 5014.96 Hz with nulls near fc ± Rs; the fine comb is the 100 Hz line spectrum of the repeating 10 ms pattern |
| `bpsk_abrupt_snr_iq_tracks.png` | 8.9 | Baseband I and Q with symbol instants marked |
| `bpsk_abrupt_snr_eye.png` | — | A clean, open eye; each capture has 15 symbols, so it is denser than the QPSK eyes |
| `bpsk_abrupt_snr_snr_sweep.png` | **8.10** | **Key plot.** Left: measured EVM (mean ± std of 50 noise draws) against SNR. Red dashed: theory including the matched-filter gain (≈125 samples/symbol ≈ 21 dB). Grey dotted: theory without that gain (reference only). Measured follows the red line about 2.5–3 dB above it, then flattens near 2.4 % at high SNR (the hardware floor). Right: BER. The measured BER is 0 at every point (shown at the plot floor) because 75 × 50 = 3750 bits per point cannot show a low BER. **Use the left panel only** |
| `ber_bpsk_curve.png` | **8.11** | **Key plot.** Simulated BER (dots, 95 % Wilson intervals) on the theory curve `Q(√(2·Eb/N0))` from 10⁻¹ to 10⁻⁸. The last three points (11–12 dB) have wide bars because they reached the 10⁸-bit limit with fewer than 100 errors |

### 7.2 Oscilloscope images ([`oscilloscope_images/`](oscilloscope_images))

Both screenshots are at the capture settings: CH1, 500 mV/div, 2.00 ms/div, 62.5 kSa/s, 2000 points.

| File | What it shows |
|---|---|
| `qpsk_nrz_scope.png` | QPSK NRZ at Point B. A constant-amplitude 1 kHz carrier with small breaks where the phase jumps. Scope readouts: Pk-Pk 1.34 V, RMS 1.889 V |
| `qpsk_rc_scope.png` | QPSK raised-cosine at Point B. The amplitude dips toward the middle at the 180° transitions (two dips visible). Pk-Pk 1.28 V, RMS 1.868 V |

The scope's **Frequency** readout (998.0 Hz and 1.053 kHz) is not a measurement of the carrier. The scope counts zero crossings, and phase jumps disturb that count. The carrier value is taken from the receiver (1000.4441 Hz).

There is no BPSK screenshot in this folder. The BPSK screenshots are in [`../spectrum/oscilloscope_images/`](../spectrum/oscilloscope_images) (`S6N.png`, `S6Nfft.png`).

---

## 8. Method details

### 8.1 Receiver chain (both scripts)

Thesis Fig. 8.1 shows this chain as a block diagram.

1. **Load.** Find the time and voltage columns in the CSV; remove the DC level.
2. **Carrier estimate.** Raise the signal to the power M (M = 2 for BPSK, 4 for QPSK). This removes the modulation and leaves a line at M·fc. Find it with a zero-padded Blackman-Harris FFT, searching only ±5 % around M·fc_nominal; divide by M. A **reliability check** warns if the peak is less than 6 dB above the median of the search band (this is what fires in 6.2).
3. **Symbol-rate correction.** The carrier and the symbol clock come from the same timer, so the symbol rate is scaled by the same clock error.
4. **Mix to baseband.** Multiply by `2·exp(−j2π·fc·t)` and low-pass filter (4th-order Butterworth, cut-off 3·Rs, zero-phase).
5. **Matched filter.** NRZ: moving average over one symbol (integrate-and-dump). Raised-cosine: none, because the transmitted pulse is already zero-ISI at the symbol centres.
6. **Symbol timing.** Oerder–Meyr estimate, then a direct search of the sampling offset that gives the lowest EVM. Symbols are read with cubic-spline interpolation; half a symbol is dropped at each end of the record.
7. **Phase and frequency clean-up.** Raise the symbols to the power M again and fit a straight line to the unwrapped phase. The slope is the leftover frequency error, the intercept the fixed rotation.
8. **Ambiguity resolution.** After step 7 the constellation can still be rotated by any multiple of 360°/M, and the capture starts at an unknown point of the repeating pattern. Every rotation and every starting position is tried against the known pattern; the one with the fewest bit errors is kept.
9. **Decision and metrics.** Nearest point (sign of I for BPSK). EVM = RMS error vector / RMS ideal amplitude; MER = −20·log₁₀(EVM). QPSK also reports I/Q gain imbalance and quadrature error.

### 8.2 Per-capture ambiguity resolution (important)

Each Single capture starts at a random point of the repeating pattern, and each capture's phase clean-up can land on a different one of the M rotations. So step 8 is done **for each capture separately**, and only the corrected symbols are pooled. The per-capture table printed by the scripts shows this: for example, the BPSK captures resolve to 180°, 0°, 180°, 180°, 0°. If step 8 is done once on all pooled symbols instead, the QPSK NRZ EVM comes out above 90 % even though every capture is fine.

### 8.3 EVM floor of the scope

The TBS1102C is 8-bit with a noise floor of about 44 dB, which limits EVM to about **0.6 %**. An EVM above this is reported as a measured value; an EVM below it (QPSK RC, 0.514 %) is reported as an upper bound, the same rule used in the SNR and THD tests.

### 8.4 EVM-vs-SNR sweep

- The noise-free basis is the **first** CSV given (`1.csv`).
- For each SNR, Gaussian noise with power `P_signal / 10^(SNR/10)` is added to the whole passband record (`P_signal` = mean square of the DC-removed capture), then the full receiver is run again. 50 noise draws per point.
- The carrier and symbol rate stay fixed (`--fc-fixed`). The sweep therefore tests the stages after carrier and timing recovery, not the carrier search itself.
- Theory: naive `EVM = 1/√SNR`, and with matched-filter gain `EVM = 1/√(SNR·sps)` with sps = 124.63 samples per symbol. The timing search uses a coarser 41-point grid inside the sweep (201 points in the main run) to keep it fast; this changes the EVM by only a few hundredths of a percent.

### 8.5 BER simulation

- Ideal receiver: perfect carrier and timing, bits → ±1, real Gaussian noise with variance `1/(2·Eb/N0)`, decision by sign.
- Stopping rule: batches of 200 000 bits until 100 errors **or** 10⁸ bits.
- 95 % Wilson score confidence interval for every point.
- Validation gate at 0 dB against `Q(√2) = 0.0786` before the sweep.
- BPSK only. QPSK (Gray-coded) has the same BER against Eb/N0, and analog PM has no bit error rate.

---

## 9. Known limitations and data notes

1. **Baseband only.** No RF up-conversion and no real radio channel.
2. **Few symbols.** QPSK: 16 symbols (32 bits) per pulse shape; BPSK: 75 bits. "0 bit errors" shows a working link; it is **not** a BER measurement. The real BER curve is the simulation.
3. **Short repeating patterns.** QPSK uses an 8-symbol pattern and BPSK a 5-bit pattern, not random data. The spectra are line spectra under the expected envelope.
4. **Incomplete QPSK pattern per capture** (Section 4), handled by pooling 3 captures.
5. **The receiver uses the known pattern** to resolve the phase ambiguity. A real receiver would use a preamble or differential coding.
6. **RC carrier taken from NRZ.** The RC captures do not give a reliable carrier estimate on their own (6.2). The NRZ value is used; both builds use the same board and clock.
7. **No repeat study.** EVM values are from 3–5 captures, without a run-to-run uncertainty budget like the frequency-accuracy test.
8. **AWGN only** in the sweep and the simulation. No fading, multipath or non-linearity.
9. **Mode 1 (two-channel baseband I/Q) not used** (Section 3.1). The PA6/DAC2 comment in `QPSK.c` does not apply to the F303RE.
10. **File names say 2.5 ms/div; the timebase was 2 ms/div** (Section 4).

---

## 10. Quick start (no hardware)

```
git clone https://github.com/moyeenurrahman/STM32_phase_modulation_research
cd STM32_phase_modulation_research/SDR
pip install numpy scipy matplotlib
mkdir reproduced_plots

# QPSK NRZ  -> EVM 1.636 %
python python/qpsk_analysis.py --csv csv/QPSK_NRZ/qpsk_mode0_nrz_2p5msdiv_cap1.csv csv/QPSK_NRZ/qpsk_mode0_nrz_2p5msdiv_cap2.csv csv/QPSK_NRZ/qpsk_mode0_nrz_2p5msdiv_cap3.csv --mode passband --shape nrz --prefix reproduced_plots/qpsk_nrz

# QPSK raised-cosine  -> EVM 0.514 % (upper bound)
python python/qpsk_analysis.py --csv csv/QPSK_raised_cosine/qpsk_mode0_rc_2p5msdiv_cap1.csv csv/QPSK_raised_cosine/qpsk_mode0_rc_2p5msdiv_cap2.csv csv/QPSK_raised_cosine/qpsk_mode0_rc_2p5msdiv_cap3.csv --mode passband --shape rc --fc 1000.4441 --fc-fixed --rs 200.0888 --prefix reproduced_plots/qpsk_rc

# NRZ vs RC figures (Fig. 8.5, 8.6, 8.7)
python python/make_nrz_rc_figures.py

# BPSK receiver + EVM-vs-SNR sweep  -> EVM 1.972 %, Fig. 8.10
python python/bpsk_analysis.py --csv csv/BPSK_from_spectral_analysis/1.csv csv/BPSK_from_spectral_analysis/2.csv csv/BPSK_from_spectral_analysis/3.csv csv/BPSK_from_spectral_analysis/4.csv csv/BPSK_from_spectral_analysis/5.csv --fc 5014.9591 --rs 501.4959 --fc-fixed --pattern 0,1,1,0,1 --shape nrz --prefix reproduced_plots/bpsk_abrupt_snr --snr-sweep --snr-min -5 --snr-max 20 --snr-step 1 --snr-realizations 50

# BPSK BER simulation  -> PASS at 0 dB, Fig. 8.11
python python/ber_simulation.py --ebno-min 0 --ebno-max 12 --ebno-step 0.5 --min-errors 100 --prefix reproduced_plots/ber_bpsk
```

The whole set runs in under a minute. Every number in Section 1 should appear in the printed output, and the files in `reproduced_plots/` should match those in `result_plots/` and `results_from_python/`.
