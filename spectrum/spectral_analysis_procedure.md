# Spectral Analysis: Measurement and Reproduction Guide

Part of the thesis *"Implementation and Performance Analysis of a Baseband Phase Modulation (PM) System with STM32 for Low-Power Wireless Communication"* (BSc EEE, RUET).

This guide explains how the **spectrum** of the STM32 baseband PM generator was measured. It tells you which firmware was flashed, how the oscilloscope was set, which CSV files were recorded, and which Python script turns each set of CSV files into a result. You can follow it to **reproduce every number and plot from the CSV files in this folder**, without any hardware.

> **Case names.** The spectral cases are called **S1 … S8**. Four of them (S1, S2, S6W, S8W) are the same captures that were made earlier for the SNR test (N-cases), so those CSV files are shared. The mapping is in Section 5.

---

## 1. What is measured, and the main results

The generator output is `s(t) = A cos( 2π fc t + kp · sin(2π fm t) )` for analog PM, or a 0°/180° phase switch for BPSK. Here `fc` = 5 kHz, `fm` = 500 Hz and the BPSK symbol rate `Rs` = 500 bit/s. The spectral analysis answers five questions:

| Question | Cases | How it is answered |
|---|---|---|
| How strong are the **DAC images** at 95 kHz and 105 kHz, and how much does the low-pass filter remove? | S1, S2, S6W | Power at 95/105 kHz relative to the carrier (dBc), before the filter (Point A) and after it (Point B) |
| Do the **PM sidebands** follow the Bessel-function pattern, and is `kp` recovered correctly? | S3, S4, S5 | Sideband amplitudes at fc ± n·fm (n = 0 … 6) are fitted to \|Jn(kp)\| |
| Does the **carrier disappear** at the Bessel null kp = 2.4048? | S2N, S4 | Carrier power in S4 compared with the unmodulated carrier in S2N |
| How wide is the signal? Does it agree with **Carson's rule**? | S3, S4, S5 | 99 % occupied bandwidth B99, compared with `2(kp + 1)·fm` |
| What does the **BPSK spectrum** look like, and does a phase ramp help? | S6N, S7N (and S6W) | Null-to-null width, first sidelobe level, harmonics |

Every case is compared with a **noise-floor capture** taken with the same scope settings (S8W for Family W, S8N for Family N), so any result that is limited by the instrument can be recognised.

### Headline results (all reproducible from this folder)

| Quantity | Result |
|---|---|
| Image at 95 / 105 kHz, **before** the filter (Point A, S1, 32 records) | **−25.69 / −26.65 dBc** (sinc theory: −25.6 / −26.5) |
| Image at 95 / 105 kHz, **after** the filter (Point B, S2, 30 records) | **−42.28 / −43.48 dBc** (an ideal filter would give −52.6 / −55.2) |
| Why the filter does less than expected | Capacitor C1 passes high frequencies around the op-amp (feedthrough). X_C1(95 kHz) = 838 Ω gives −16.5 dB; the measured filter contribution is **−16.59 dB**. The MCP6002 gain-bandwidth product is only 1 MHz |
| Recovered `kp` (S3, S4, S5) | **1.0000, 2.4020, 3.8280** rad for commanded 1.0000, 2.4048, 3.8317 (errors 0.00 %, −0.12 %, −0.10 %) |
| Bessel-pattern fit residual | 0.0021, 0.0014, 0.0016 |
| Carrier null at kp = 2.4048 (S4 against S2N) | **54.90 dB** below the unmodulated carrier |
| Sidebands at kp = 3.8317 (S5) | First sidebands at −45.8 dBc; the second sidebands equal the carrier (**−0.00 dBc**) |
| B99 (S3 / S4 / S5), band 2 – 8 kHz | **2031 / 3094 / 5031 Hz** against Carson 2000 / 3405 / 4832 Hz (**+1.6 / −9.1 / +4.1 %**) |
| B99 of the unmodulated carrier (S2N) | 94 Hz |
| BPSK null-to-null width, abrupt (S6N) | **1031 Hz** (ideal 1000 Hz) |
| BPSK null-to-null width, phase ramp (S7N) | **1281 Hz**. The ramp made the spectrum **wider**, not narrower |
| BPSK first sidelobe | −11.60 dBc (S6N), −11.37 dBc (S7N), against −13.3 dBc for ideal rectangular NRZ |
| BPSK 2nd harmonic | −34.54 dBc (S6N), −29.07 dBc (S7N): the ramp is 5.5 dB **worse** |

---

## 2. Equipment and wiring

| Item | Use |
|---|---|
| STM32 Nucleo-F303RE | Waveform generator: DAC1_OUT1 on **PA4**, TIM6 trigger → circular DMA |
| Sallen-Key Butterworth low-pass filter (MCP6002 at 4.05 V, R = 5.6 kΩ, C1 = 2 nF, C2 = 1 nF, f_c ≈ 20.1 kHz) | Reconstruction filter |
| Tektronix TBS1102C oscilloscope (TBS1000C series) | CH1 only, 8-bit, **fixed 2000-point record** |
| ×1 probe | All spectral captures |
| USB flash drive | CSV export |
| STM32CubeIDE | Build and flash |
| Python 3 | numpy, scipy, matplotlib |

**Two test points:**

- **Point A** = raw DAC output (PA4), before the filter. Shows the staircase and the images.
- **Point B** = filter output. Used for everything except the "before" half of the image comparison.

Use one common ground and keep the probe ground clip short.

---

## 3. Firmware (STM32 code)

Two firmware files were used. Both are the `main.c` of a CubeMX project and generate the same signals. They differ only in the test modes they offer.

| File | Used for |
|---|---|
| [`firmware/spectral_analysis.c`](firmware/spectral_analysis.c) | **Family N**: cases S2N, S3, S4, S5, S6N, S7N, S8N |
| [`../SNR/firmware/SNR.c`](../SNR/firmware/SNR.c) | **Family W**: cases S1, S2, S6W, S8W. These captures were made during the SNR test, so the SNR firmware generated them |

### 3.1 What the firmware does

- Clock: **HSI 8 MHz, no PLL**. The timer reload value is computed at run time. Use the same clock as in the SNR and frequency-accuracy tests.
- Fills a **1000-sample buffer** with `2048 + 1200 · cos(...)` and streams it to the DAC with **TIM6 trigger → circular DMA**. The CPU does nothing while the signal runs.
- DAC update rate **100 kSa/s**. The buffer holds 50 carrier cycles, 5 cycles of the 500 Hz message, or 5 BPSK symbols, so it repeats every 10 ms (**100 Hz**). The output swing is about 0.68 – 2.62 V.
- Because the buffer repeats, the BPSK signal is a **periodic 5-bit pattern** (`0x16` = 0, 1, 1, 0, 1). Its spectrum is a line spectrum with 100 Hz spacing, not the continuous spectrum of random data (see Section 9).
- At start-up the LED flashes **SPEC_MODE + 1** times, so you can check which build is running.

### 3.2 Choosing the test: one line to change

In `spectral_analysis.c`, change this line, then **Clean → Build → Flash**:

```c
#define SPEC_MODE   MODE_CARRIER
```

| `SPEC_MODE` | LED flashes | Case | Signal |
|---|---|---|---|
| `MODE_FLOOR` | 1 | **S8N** | DMA runs, buffer is all mid-code (2048). Instrument noise floor |
| `MODE_CARRIER` | 2 | **S2N** | Unmodulated 5 kHz carrier |
| `MODE_PM_1` | 3 | **S3** | Analog PM, kp = 1.0 |
| `MODE_PM_BESSEL` | 4 | **S4** | Analog PM, kp = 2.4048 (first zero of J0, carrier nulls) |
| `MODE_PM_NULL1` | 5 | **S5** | Analog PM, kp = 3.8317 (first zero of J1, first sidebands null) |
| `MODE_BPSK_ABRUPT` | 6 | **S6N** | BPSK, hard 0°/180° switching |
| `MODE_BPSK_RAMP` | 7 | **S7N** | BPSK, phase ramp: `RAMP_SAMPLES 60`, `RAMP_STEPS 6` (the phase moves in 30° steps over the first 0.6 ms of each 2 ms symbol) |
| `MODE_DC_STATIC` | 8 | (optional) | DMA stopped, DAC held at 2048. Not used in the results |

For Family W, the SNR firmware modes were used (`SNR_MODE` in `SNR.c`):

| `SNR_MODE` | Case | SNR name |
|---|---|---|
| `MODE_DC_DMA` | **S8W** (noise floor) | N1 |
| `MODE_CARRIER` | **S1** (Point A) and **S2** (Point B) | N2 and N3 |
| `MODE_BPSK` | **S6W** | N6 |

---

## 4. Oscilloscope settings

The TBS1000C always saves **2000 points**. So the timebase decides both the sample rate and the frequency resolution. Two settings are used, called **Family W** (wide) and **Family N** (narrow).

| | Family W | Family N |
|---|---|---|
| Timebase | **400 µs/div** | **2 ms/div** |
| Sample rate | 312.5 kSa/s | 62.5 kSa/s |
| Sample interval (in the CSV header) | 3.2 µs | 16 µs |
| Record length | 6.4 ms | 32 ms |
| FFT bin width | 156.25 Hz | 31.25 Hz |
| Highest frequency (Nyquist) | 156.25 kHz | 31.25 kHz |
| Used for | DAC images at 95/105 kHz, before and after the filter | Sidebands, Bessel fit, B99, harmonics, BPSK lobes |
| Test points | A and B | **B only** |

**Why Family N is Point B only.** At 62.5 kSa/s the 105 kHz image folds to **20 kHz** (the 4th harmonic) and the 95 kHz image folds to **30 kHz**. At Point B these images are small and harmless. At Point A they would sit on top of the harmonics. The script refuses to analyse a Point A capture in Family N.

**Settings common to all captures**

| Setting | Value |
|---|---|
| Vertical | **500 mV/div**, ×1 probe, DC coupling |
| Acquisition | **Hi-Res** mode |
| Export | USB flash drive. (The OpenChoice PC software cannot read Hi-Res captures, so the pendrive was used) |

**Reading the CSV header.**

```
Record Length,2000,Points
Sample Interval,16.0000000E-6,s       <- the number the scripts use
Horizontal Scale,3.19999999E-3,       <- NOT the front-panel timebase
```

The `Horizontal Scale` in the header is the whole 2000-point record divided by 10. It is **not** the timebase knob setting. (Family N shows 3.2 ms in the header for 2 ms/div; Family W shows 640 µs for 400 µs/div.) Use `Sample Interval`. The script computes the sample rate from the time column and rejects any file that is more than 1 % away from the expected rate.

---

## 5. Folder layout, file names and case table

```
spectrum/
├── csv/                       oscilloscope CSV files, one folder per case
│   ├── S1/ S2/ S6W/ S8W/                  Family W (400 µs/div)
│   └── S2N/ S3/ S4/ S5/ S6N/ S7N/ S8N/    Family N (2 ms/div)
├── firmware/spectral_analysis.c
├── python/
│   ├── analyse_spectrum.py          main analysis, all cases
│   ├── figure_image_theory.py       -> image_theory.png
│   ├── figure_bessel_fit.py         -> bessel_fit.png
│   ├── figure_b99_carson.py         -> b99_carson.png
│   └── figure_bpsk_comparison.py    -> bpsk_comparison.png
├── oscilloscope_images/       screenshots from the scope (Section 7.2)
├── result_plots/              plots used in the thesis (Section 7.1)
├── reproduced_plots/          created when you run the scripts
└── spectral_analysis_procedure.md   this guide
```

**File-name rule.** Inside each case folder the files are simply `1.csv`, `2.csv`, … one file per capture. The folder name is the case name.

| Case | Family | Point | Firmware mode | CSV files | Same files as |
|---|---|---|---|---|---|
| **S1** | W | A | `SNR.c` `MODE_CARRIER` | 32 | `SNR/csv/N2` |
| **S2** | W | B | `SNR.c` `MODE_CARRIER` | 32 (2 are duplicates, so 30 are used) | `SNR/csv/N3` |
| **S6W** | W | B | `SNR.c` `MODE_BPSK` | 8 | `SNR/csv/N6` |
| **S8W** | W | B | `SNR.c` `MODE_DC_DMA` | 32 | `SNR/csv/N1` |
| **S2N** | N | B | `MODE_CARRIER` | 5 | new capture |
| **S3** | N | B | `MODE_PM_1` | 5 | new capture |
| **S4** | N | B | `MODE_PM_BESSEL` | 5 | new capture |
| **S5** | N | B | `MODE_PM_NULL1` | 5 | new capture |
| **S6N** | N | B | `MODE_BPSK_ABRUPT` | 5 | new capture |
| **S7N** | N | B | `MODE_BPSK_RAMP` | 5 | new capture |
| **S8N** | N | B | `MODE_FLOOR` | 5 | new capture |

"Same files as" was checked by comparing file hashes (MD5). Family W was not captured again for the spectral analysis. The SNR captures were reused, so the SNR chapter and the spectral chapter report the same image values.

In S2, `2.csv` and `3.csv` are byte-identical, and so are `4.csv` and `6.csv` (a capture was saved twice). The script detects this by MD5 and counts each file only once.

---

## 6. Reproducing the results

### 6.0 Setup

```
pip install numpy scipy matplotlib
```

Tested with numpy 2.5, scipy 1.18 and matplotlib 3.11. No staging or renaming is needed. The scripts read `spectrum/csv/<case>/` directly and write all plots to `spectrum/reproduced_plots/`, so they can be run from any folder.

On a computer without a screen, set `MPLBACKEND=Agg` first (Linux/macOS: `export MPLBACKEND=Agg`; Windows: `set MPLBACKEND=Agg`). Otherwise the four figure scripts open a window at the end, which you close to finish.

### 6.1 Main analysis

```
cd spectrum/python
python analyse_spectrum.py
```

The script does these steps for every case:

1. **Screen the files.** MD5 duplicates are dropped. Files whose sample rate is more than 1 % from the family value, or that have fewer than 1500 points, are rejected.
2. **Remove the DC level**, then compute one **periodogram over the full 2000-point record** with a **Blackman-Harris** window.
3. **Average the records in the power domain** (not the amplitude domain).
4. Run the analysis for the case type (Section 8) and print the numbers.
5. Save one PSD plot per case (`<case>_psd.png/.pdf`) and the before/after filter plot.

**Expected file-screening table:**

```
Case   Family  Point  Found   Dupes   Bad    Used   Bin (Hz)
S1     W       A      32      0       0      32     156.25
S2     W       B      32      2       0      30     156.25
S6W    W       B      8       0       0      8      156.25
S8W    W       B      32      0       0      32     156.25
S2N    N       B      5       0       0      5      31.25
S3     N       B      5       0       0      5      31.25
S4     N       B      5       0       0      5      31.25
S5     N       B      5       0       0      5      31.25
S6N    N       B      5       0       0      5      31.25
S7N    N       B      5       0       0      5      31.25
S8N    N       B      5       0       0      5      31.25
```

Report the **Used** count in a write-up, not the number of files.

**Key numbers the script must print:**

| Case | Printed result |
|---|---|
| S1 | Image 95 kHz −25.69 dBc, 105 kHz −26.65 dBc |
| S2 | Image 95 kHz −42.28 dBc (+20.3 dB above the floor), 105 kHz −43.48 dBc (+18.3 dB above the floor) |
| S6W | Image 95 kHz −43.00 dBc, 105 kHz −44.06 dBc |
| S2N | Carrier −3.30 dB(V²); B99 = 94 Hz; harmonics H2 … H6 = −47.84, −54.82, −47.68, −54.89, −58.98 dBc |
| S3 | Recovered kp 1.0000 (residual 0.0021); B99 2031 Hz (Carson 2000, +1.6 %) |
| S4 | Recovered kp 2.4020 (−0.12 %, residual 0.0014); **Carrier null depth vs S2N carrier: 54.90 dB**; B99 3094 Hz (Carson 3405, −9.1 %) |
| S5 | Recovered kp 3.8280 (−0.10 %, residual 0.0016); n = 2 sidebands −0.00 dBc; B99 5031 Hz (Carson 4832, +4.1 %) |
| S6N | Lower/upper null 4500 / 5531 Hz, so **1031 Hz**; upper/lower first sidelobe −11.60 / −10.44 dBc; H2 −34.54 dBc |
| S7N | Lower/upper null 4375 / 5656 Hz, so **1281 Hz**; upper/lower first sidelobe −11.37 / −11.55 dBc; H2 −29.07 dBc |

### 6.2 Figure scripts

Each of the four figure scripts imports `analyse_spectrum.py`, so run them from the `python/` folder. They recompute their numbers from the CSV files and print them.

```
python figure_image_theory.py        # image_theory.png
python figure_bessel_fit.py          # bessel_fit.png
python figure_b99_carson.py          # b99_carson.png
python figure_bpsk_comparison.py     # bpsk_comparison.png
```

| Script | Prints |
|---|---|
| `figure_image_theory.py` | `Used 30 averaged records for S2 (found 32, 2 duplicates)`; ideal prediction −52.6 / −55.2 dB; feedthrough model −42.1 / −43.9 dB |
| `figure_bessel_fit.py` | S3 1.0000 / 1.0000 / 0.0021, S4 2.4048 / 2.4020 / 0.0014, S5 3.8317 / 3.8280 / 0.0016 (commanded / recovered / residual) |
| `figure_b99_carson.py` | S3 2031 vs 2000 Hz (+1.6 %), S4 3094 vs 3405 Hz (−9.1 %), S5 5031 vs 4832 Hz (+4.1 %), all PASS (±20 %) |
| `figure_bpsk_comparison.py` | Null-to-null 1031 Hz (S6N) and 1281 Hz (S7N); H2 −34.54 / −29.07 dBc |

The plot titles say "Figure 7.x" because those are the figure numbers in the thesis.

### 6.3 Output files

All outputs go to `spectrum/reproduced_plots/`. The images in `result_plots/` were made from the same data. The numbers are identical; small differences in fonts or spacing can appear between matplotlib versions.

| Script | Writes (PNG and PDF) | Same as in `result_plots/` |
|---|---|---|
| `analyse_spectrum.py` | `S1_psd`, `S2_psd`, `S2N_psd`, `S3_psd`, `S4_psd`, `S5_psd`, `S6N_psd`, `S6W_psd`, `S7N_psd` | same names |
| `analyse_spectrum.py` | `before_after_LPF` | `before_after_LPF.png` |
| `figure_image_theory.py` | `image_theory` | `image_theory.png` |
| `figure_bessel_fit.py` | `bessel_fit` | `bessel_fit.png` |
| `figure_b99_carson.py` | `b99_carson` | `b99_carson.png` |
| `figure_bpsk_comparison.py` | `bpsk_comparison` | `bpsk_comparison.png` |

The two floor cases (S8W, S8N) have no plot of their own. They appear as the thin orange "instrument floor" line in the other plots.

---

## 7. What each plot and screenshot shows

### 7.1 Result plots ([`result_plots/`](result_plots))

"Fig." is the figure number in the thesis (Section 7.8). All PSD plots show power spectral density in dB(V²/Hz) against frequency, averaged over the records in the "Used" column. In the single-case plots the orange line is the instrument noise floor of the same family. Anything close to this line is limited by the scope, not by the generator.

**Plots made by their own scripts**

| File | Fig. | CSV used | What it shows and what to look for |
|---|---|---|---|
| `before_after_LPF.png` | 7.5 | S1, S2, S8W | **Key plot.** Spectrum before (blue) and after (orange) the filter, 0 – 156 kHz, with the noise floor (green dotted). The 5 kHz carrier is the same in both. The strongest peaks near 95 and 105 kHz (dashed lines) are the DAC images: they are tall before the filter and about 17 dB lower after it, but not removed |
| `image_theory.png` | 7.6 | S2 | Measured spectrum at Point B relative to the carrier (blue) against the **ideal** curve (sinc × 2nd-order Butterworth, green dashed). The ideal curve has a deep notch at 100 kHz. The red crosses are the feedthrough estimate at 95/105 kHz (−42.1 / −43.9 dB), which matches the measured images; the ideal curve would put them 10 dB lower |
| `bessel_fit.png` | 7.10 | S3, S4, S5 | **Key plot.** Three panels, one per kp. Bars: measured sideband amplitude (blue) and the fitted Bessel prediction (orange) for orders n = 0 … 6, in dB. The panel titles give commanded kp, recovered kp and the fit residual. At kp = 2.4048 the n = 0 bar (carrier) is very low; at kp = 3.8317 the n = 1 bar is very low |
| `b99_carson.png` | 7.11 | S3, S4, S5 | Measured B99 (blue) against Carson's rule (orange). Error bars are ±20 %. The percentage above each pair is the difference. All three fall inside the band |
| `bpsk_comparison.png` | 7.14 | S6N, S7N | Left: the BPSK main lobe around 5 kHz, abrupt (blue) against ramped (orange); dotted lines mark the measured nulls (the ramped nulls are further apart). Right: 0 – 31 kHz with harmonics H2 … H6 marked |

**Per-case PSD plots** (all made by `analyse_spectrum.py`)

| File | Fig. | CSV used | What it shows |
|---|---|---|---|
| `S1_psd.png` | 7.3 | S1 | Point A, 0 – 156 kHz, 95/105 kHz marked. The raw DAC output, before the filter: strong images |
| `S2_psd.png` | 7.4 | S2 | Point B, 0 – 156 kHz. The same signal after the filter: weaker images |
| `S2N_psd.png` | (not numbered) | S2N | Unmodulated 5 kHz carrier, 0 – 25 kHz, B99 band shaded. A single narrow line (B99 = 94 Hz) plus the harmonics. This is the reference for the carrier null in S4 |
| `S3_psd.png` | 7.7 | S3 | kp = 1: the carrier and sideband pairs at 500 Hz spacing, falling quickly with order (about −4.8, −16.5, −31.9, −48 dBc for n = 1 … 4). B99 band shaded |
| `S4_psd.png` | 7.8 | S4 | kp = 2.4048: **the carrier at 5 kHz is almost gone** (54.90 dB below S2N) while the sidebands stay strong |
| `S5_psd.png` | 7.9 | S5 | kp = 3.8317: the **first sidebands at 4.5 / 5.5 kHz are almost gone**; the carrier and the second sidebands are equal in height. The widest signal (B99 5031 Hz) |
| `S6N_psd.png` | 7.12 | S6N | BPSK, abrupt switching. A main lobe with nulls near 4.5 and 5.5 kHz. The fine comb is the 100 Hz line spectrum of the repeating 5-bit pattern |
| `S7N_psd.png` | 7.13 | S7N | BPSK with the phase ramp. A wider main lobe (nulls at 4375 and 5656 Hz) |
| `S6W_psd.png` | (not numbered) | S6W | BPSK at Point B on the wide span, 0 – 156 kHz. The broad lobe is centred at 5 kHz and falls off with frequency. The images near 95 and 105 kHz are at −43.00 / −44.06 dBc |

The PSD is on a log scale. A flat band near −85 to −90 dB(V²/Hz) (Family N) is the 8-bit scope floor, not signal.

### 7.2 Oscilloscope images ([`oscilloscope_images/`](oscilloscope_images))

These are **screen captures** from the scope itself. They show the setup and give a quick visual check. **The numbers in this guide come from the CSV files, not from these screenshots.**

| File | Case | What it shows |
|---|---|---|
| `S2N.png` | S2N | Time-domain view of the unmodulated 5 kHz carrier at 2 ms/div, 62.5 kSa/s. Peak-to-peak 1.96 V, maximum 2.63 V |
| `S3.png` | S3 | The scope's own FFT of the kp = 1 signal: carrier and sideband pairs |
| `S4.bmp` | S4 | Time-domain view at kp = 2.4048, 1 ms/div, 125 kSa/s. (BMP format; any image viewer opens it) |
| `S4fft.png` | S4 | Scope FFT at kp = 2.4048: the centre (carrier) peak is replaced by a dip between the sideband peaks |
| `S5.png` | S5 | Time-domain view at kp = 3.8317 |
| `S5fft.png` | S5 | Scope FFT at kp = 3.8317 |
| `S6N.png` | S6N | Time-domain view of BPSK with abrupt switching: the carrier phase reverses at the symbol changes |
| `S6Nfft.png` | S6N | Scope FFT of the abrupt BPSK |
| `S7N.png` | S7N | Time-domain view of BPSK with the phase ramp |
| `S7Nfft.png` | S7N | Scope FFT of the ramped BPSK |

Notes on these images:

- **The FFT screenshots are qualitative.** `S3.png`, `S4fft.png`, `S5fft.png`, `S6Nfft.png` and `S7Nfft.png` were taken at 400 µs/div (about 312 kSa/s) with the scope's built-in **Hanning** window and 20 dB/div. The CSV files for these cases (Family N) were taken at 2 ms/div. So the screenshots cannot be compared bin by bin with the analysis plots. Use them to see the carrier null (S4) and the sideband pattern on the instrument itself.
- **Instrument clock.** The date drawn on the scope screen (10 Sep 2020) and the date in the screenshot caption strip and CSV header note (for example `15-Aug-26`) do not agree. Do not use the on-screen date as the capture date.
- **Missing screenshots.** None exist for S1, S2, S6W, S8N and S8W in this folder. The equivalent images for the shared Family W captures are in `SNR/oscilloscope_images/`:

| SNR image | Case | Content |
|---|---|---|
| `N2N3.png` | S1 and S2 | Time domain: CH1 filtered sine, CH2 raw DAC staircase |
| `N2fft.png` | S1 | Scope FFT at Point A, with the image peaks |
| `N3.png`, `N3fft.png` | S2 | Point B time domain and scope FFT |
| `N6_scope.png` | S6W | BPSK capture |

---

## 8. Method details

- **Spectrum estimate.** The mean is removed and one periodogram is computed over all 2000 points (`scipy.signal.welch` with `nperseg` equal to the record length and no overlap) using a Blackman-Harris window. A single periodogram is used because a 2000-point record cannot be split into segments without losing the 31.25 Hz resolution. Blackman-Harris is used because the scope clock and the STM32 clock are independent, so the signal never sits exactly on a bin. Records are averaged in the **power** domain.
- **Line power.** The power of a line is the **sum of the bins within ±3 bins** of the expected frequency, not just the peak bin, so the answer does not depend on where the line falls between bins.
- **Image suppression (Family W).** Power at 95 and 105 kHz relative to the carrier, in dBc. If an image is within 3 dB of the instrument floor, the result is reported only as a limit ("at least"). All images here are 17 dB or more above the floor, so they are real measurements.
- **Image prediction.** Ideal: the DAC zero-order-hold `sinc` response times the 2nd-order Butterworth response. The measured Point B images are higher because of C1 feedthrough (X_C1 at 95 kHz = 838 Ω), which the ideal model does not include.
- **Bessel recovery of kp (Family N).** The amplitudes at fc + n·fm and fc − n·fm are averaged for n = 1 … 6, and the carrier gives n = 0. This 7-point pattern is normalised and fitted by least squares to |Jn(kp)| for n = 0 … 6. The fit residual shows how Bessel-like the pattern is. This is used instead of the J1/J0 ratio, which diverges at kp = 2.4048 (J0 = 0).
- **Occupied bandwidth.** The PSD is integrated over a **symmetric 2 – 8 kHz** band around the 5 kHz carrier, and B99 is the width that contains 99 % of that power. The band is symmetric and stops below the 2nd harmonic at 10 kHz. Carson's rule for PM is `B = 2(kp + 1)·fm`.
- **BPSK.** B99 is **not** reported (the true B99 of rectangular NRZ is about 20·Rs = 10 kHz, outside the band). Instead the nulls and the first sidelobe near fc ± Rs are found, and the harmonics H2 … H6 are measured. No SINAD or ENOB is reported for BPSK.
- **Floor check.** The median PSD of each case is compared with the floor capture of the same family. A margin below 3 dB means the result is instrument-limited.

---

## 9. Known limitations and data notes

1. **Periodic BPSK pattern.** The firmware repeats a 5-bit pattern every 10 ms, so the BPSK spectrum is a line spectrum (100 Hz spacing) under a sinc-like envelope. The null-to-null width and sidelobe level follow the envelope, but the lines are not the spectrum of random data.
2. **8-bit scope.** The floor of the Family N captures is about −87 dB(V²/Hz). Weak lines at the level of the S3 n = 5 and 6 sidebands (about −55 dBc) are close to it, so the Bessel fit uses the whole 7-point pattern and not single ratios.
3. **S4 reference values.** The harmonic and shoulder values printed for S4 in dBc are referenced to the **nulled carrier**. They are large positive numbers and are **not meaningful**. For S4 use the null depth (54.90 dB against S2N) and the Bessel fit. The script prints a note about this.
4. **Frequency resolution.** The bin width is 31.25 Hz (Family N), so mains hum at 50 Hz cannot be separated from the carrier skirts.
5. **Unidentified spurs.** The unmodulated carrier (S2N) shows small lines near 3.0, 3.5, 6.5 and 8.0 kHz at about −50 to −60 dBc. Their source was not identified. They fall on sideband positions of the PM cases. At n = 3 (3.5 / 6.5 kHz) the real S3 sidebands (−31.9 dBc) are about 17 dB or more above them. At n = 4 and beyond the S3 sidebands (−48 dBc and lower) are at a similar level, so those amplitudes carry some uncertainty. The recovered kp stays within 0.12 % because the strong low-order sidebands dominate the fit.
6. **Duplicates in S2.** Two files were saved twice. After removal 30 records are averaged, against 32 for S1. This changes the variance of the estimate, not the mean.
7. **S4 B99 is 9.1 % below Carson's rule.** This is inside the ±20 % tolerance used here. Carson's rule is an approximation, not an exact value.
8. **The phase ramp did not help.** The ramp (30° steps over 0.6 ms) creates a short frequency offset during each transition. That widens the main lobe (1281 Hz against 1031 Hz) and raises H2 by 5.5 dB. This is reported as it is: this ramp design does not reduce the spectral spread.
9. **Screenshots and CSV differ in timebase** (see Section 7.2). Only the CSV files are used for numbers.
10. **Family W data were not captured again.** S1, S2, S6W and S8W were recorded with the SNR firmware (`SNR.c`), not with `spectral_analysis.c`.
11. **One operating point.** One filter, one carrier (5 kHz) and one DAC amplitude (1200 codes). No temperature or supply-voltage sweep.

---

## 10. Quick start (no hardware)

```
git clone https://github.com/moyeenurrahman/STM32_phase_modulation_research
cd STM32_phase_modulation_research/spectrum/python
pip install numpy scipy matplotlib

python analyse_spectrum.py          # all cases: screening table, numbers, PSD plots
python figure_image_theory.py       # Fig 7.6
python figure_bessel_fit.py         # Fig 7.10
python figure_b99_carson.py         # Fig 7.11
python figure_bpsk_comparison.py    # Fig 7.14
```

The plots appear in `spectrum/reproduced_plots/`. Every number in Section 1 should appear in the printed output.
