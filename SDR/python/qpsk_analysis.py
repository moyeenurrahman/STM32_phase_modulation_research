#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
qpsk_analysis.py
================

QPSK constellation recovery and EVM analysis for STM32 baseband captures.

Project : Implementation and Performance Analysis of a Baseband Phase
          Modulation System with STM32 for Low-Power Wireless Communication
Board   : NUCLEO-F303RE, DAC1 + TIM6 TRGO + circular DMA
Scope   : Tektronix TBS1000C, 2000-point record, 8-bit vertical

Companion to the QPSK firmware (main.c).  Handles both firmware modes:

  MODE 0 (passband)  One scope channel on PA4.  The script recovers I and Q
                     by mixing the capture down with a locally generated
                     carrier, exactly as a software receiver would.

  MODE 1 (baseband)  Two scope channels, I on PA4 and Q on PA5.  The script
                     reads them directly.

Processing chain (passband mode)
--------------------------------
  1. Load and de-trend the capture
  2. Estimate the true carrier frequency by the 4th-power method
  3. Scale the symbol rate by the same clock error
  4. Mix down to complex baseband, low-pass filter
  5. Matched filter (integrate-and-dump for NRZ)
  6. Symbol timing recovery: Oerder-Meyr, refined by minimum-EVM search
  7. Residual frequency and phase correction from the symbol phases
  8. Resolve the 4-fold phase ambiguity against the known bit pattern
  9. Compute EVM, MER, I/Q imbalance, and bit errors

Typical use
-----------
  # Passband, single capture
  python qpsk_analysis.py --csv capture_2p5ms.csv --mode passband

  # Passband, several captures pooled for better EVM statistics
  python qpsk_analysis.py --csv cap1.csv cap2.csv cap3.csv --mode passband

  # Baseband I/Q, both channels in one file
  python qpsk_analysis.py --csv iq_capture.csv --mode baseband

  # Baseband I/Q, one file per channel
  python qpsk_analysis.py --csv chI.csv --csv-q chQ.csv --mode baseband

Recommended scope settings
--------------------------
  2.5 ms/div  ->  50.0 kSa/s, 40 ms record,  8 symbols  (one full pattern)
  5.0 ms/div  ->  25.0 kSa/s, 80 ms record, 16 symbols  (better statistics)

  The TBS1000C rule is fs = 125 / (timebase in ms) kSa/s with a fixed
  2000-point record.  Avoid going slower than 5 ms/div: the residual DAC
  image near 99 kHz can alias back into the measurement band once the
  scope sample rate drops far below it.

Instrument floor
----------------
  The TBS1000C digitiser is 8-bit.  The measured noise floor of about
  44 dB sets an EVM floor near 0.6 percent.  Any EVM below roughly 1
  percent should be reported as instrument-limited, not as a property of
  the STM32 transmitter.

Requires: numpy, scipy, matplotlib
"""

import argparse
import csv
import os
import sys

import numpy as np
from scipy import signal
from scipy.interpolate import CubicSpline

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ===========================================================================
#  DEFAULTS - these match the firmware CONFIG block
# ===========================================================================

FC_NOMINAL_HZ = 1000.0    # carrier frequency set in firmware
RS_NOMINAL_HZ = 200.0     # symbol rate set in firmware

# Gray-coded QPSK map used by the firmware.
#   d = 0 (00) -> I = +1, Q = +1   ( 45 deg)
#   d = 1 (01) -> I = -1, Q = +1   (135 deg)
#   d = 2 (10) -> I = +1, Q = -1   (315 deg)
#   d = 3 (11) -> I = -1, Q = -1   (225 deg)
QPSK_I = np.array([+1.0, -1.0, +1.0, -1.0])
QPSK_Q = np.array([+1.0, +1.0, -1.0, -1.0])

# Symbol pattern programmed into the firmware
TX_PATTERN = np.array([0, 1, 3, 2, 0, 3, 1, 2], dtype=int)

# Ideal constellation, unit average power
IDEAL = (QPSK_I + 1j * QPSK_Q) / np.sqrt(2.0)


# ===========================================================================
#  CSV LOADING
# ===========================================================================

def _to_float(s):
    """Parse a CSV field as a float, or return NaN."""
    try:
        return float(s.strip())
    except (ValueError, AttributeError):
        return np.nan


def load_tek_csv(path, verbose=True):
    """
    Load a Tektronix CSV export and return (t, [ch1, ch2, ...]).

    Written to survive the several formats these scopes produce:

      * OpenChoice Desktop 'Get Data' saves, which put metadata such as
        Record Length and Sample Interval in the first two columns and the
        waveform in the right-hand columns.
      * Direct USB pendrive saves.
      * Plain TIME,CH1[,CH2] files.

    The method is format-agnostic.  It looks for the longest run of rows
    whose right-hand fields all parse as numbers, then identifies the time
    axis as the column with a near-constant positive step.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(path)

    with open(path, "r", newline="", encoding="utf-8", errors="replace") as f:
        rows = [r for r in csv.reader(f)]

    if not rows:
        raise ValueError(f"{path}: file is empty")

    ncols = max(len(r) for r in rows)

    # Build a float matrix, padding short rows with NaN
    mat = np.full((len(rows), ncols), np.nan)
    for i, r in enumerate(rows):
        for j, cell in enumerate(r):
            mat[i, j] = _to_float(cell)

    # Find the longest block of consecutive rows that share the same
    # pattern of numeric columns.  That block is the waveform data.
    numeric_mask = ~np.isnan(mat)
    signatures = [tuple(row) for row in numeric_mask]

    best_start, best_len = 0, 0
    i = 0
    while i < len(signatures):
        j = i
        while j + 1 < len(signatures) and signatures[j + 1] == signatures[i]:
            j += 1
        run = j - i + 1
        # Require at least two numeric columns, else it is a metadata block
        if sum(signatures[i]) >= 2 and run > best_len:
            best_start, best_len = i, run
        i = j + 1

    if best_len < 10:
        raise ValueError(f"{path}: could not find a waveform data block")

    block = mat[best_start:best_start + best_len, :]
    cols = np.where(~np.isnan(block[0]))[0]
    data = block[:, cols]

    # Identify the time column: monotonic, with a near-constant step.
    #
    # Tolerance note: Tektronix OpenChoice CSV exports write timestamps as
    # text with ~9 significant figures (e.g. "-15.6160006E-3"). Rounding at
    # that precision alone produces a step-to-step relative variation on the
    # order of 1e-5, which is a property of the export format, not evidence
    # that the column isn't a time axis. A 1e-6 relative tolerance rejects
    # every real capture; 1e-2 (1 percent) is loose enough to accept normal
    # export rounding while still rejecting an oscillating signal column,
    # which fails the strict monotonicity check above regardless of this
    # tolerance.
    time_col = None
    for k in range(data.shape[1]):
        col = data[:, k]
        d = np.diff(col)
        if np.all(d > 0) and np.std(d) < 1e-2 * abs(np.mean(d)) + 1e-15:
            time_col = k
            break

    if time_col is None:
        # No time axis found. Fall back to sample index and look for a
        # Sample Interval entry in the metadata.
        dt = None
        for r in rows[:best_start]:
            joined = ",".join(r).lower()
            if "sample interval" in joined:
                for cell in r:
                    v = _to_float(cell)
                    if not np.isnan(v) and v > 0:
                        dt = v
                        break
        if dt is None:
            raise ValueError(
                f"{path}: no time column and no Sample Interval in header. "
                "Re-export with the time axis, or pass --fs explicitly."
            )
        t = np.arange(data.shape[0]) * dt
        channels = [data[:, k] for k in range(data.shape[1])]
    else:
        t = data[:, time_col]
        channels = [data[:, k] for k in range(data.shape[1]) if k != time_col]

    if verbose:
        fs = 1.0 / np.mean(np.diff(t))
        print(f"  {os.path.basename(path)}: {len(t)} points, "
              f"fs = {fs / 1e3:.3f} kSa/s, "
              f"span = {(t[-1] - t[0]) * 1e3:.2f} ms, "
              f"{len(channels)} channel(s)")

    return t, channels


# ===========================================================================
#  CARRIER ESTIMATION
# ===========================================================================

def estimate_carrier(x, fs, fc_nominal, search_frac=0.05):
    """
    Estimate the true carrier frequency by the 4th-power method.

    A QPSK signal carries no discrete line at fc, because the four phases
    average to zero.  Raising the signal to the fourth power collapses all
    four phase states onto one, because 4 x 45 deg, 4 x 135 deg, 4 x 225 deg
    and 4 x 315 deg are all equal to 180 deg modulo 360.  That leaves a
    clean tone at 4 x fc, which can be located directly.

    The search is restricted to a window around 4 x fc_nominal so that DAC
    images and filter artefacts cannot be mistaken for the carrier.

    Returns the estimated fc in Hz.
    """
    n = len(x)

    # Window before the nonlinearity, or the record edges dominate the FFT.
    w = signal.windows.blackmanharris(n)
    x4 = (x * w) ** 4

    # Heavy zero-padding buys interpolation resolution, not information.
    nfft = int(2 ** np.ceil(np.log2(n * 16)))
    X = np.fft.rfft(x4 - np.mean(x4), n=nfft)
    freqs = np.fft.rfftfreq(nfft, d=1.0 / fs)
    mag = np.abs(X)

    f_target = 4.0 * fc_nominal
    lo = f_target * (1.0 - search_frac)
    hi = f_target * (1.0 + search_frac)

    band = (freqs >= lo) & (freqs <= hi)
    if not np.any(band):
        raise ValueError(
            f"4*fc = {f_target:.1f} Hz falls outside the captured bandwidth "
            f"(fs = {fs:.0f} Sa/s). Use a faster timebase."
        )

    idx = np.argmax(np.where(band, mag, 0.0))

    # Warn if the located peak has poor prominence relative to the rest of
    # the search band, rather than checking literal bin position. A window
    # like Blackman-Harris has a wide main lobe, so when the true 4*fc line
    # is weak or absent, leaked energy can produce a local maximum several
    # bins inside the window edge, not necessarily at the edge itself --
    # checking array-index proximity to the boundary misses that case.
    # Prominence -- how far the chosen peak stands above the typical level
    # across the rest of the search band -- is a more direct test of
    # whether a genuine spectral line was actually found. This is exactly
    # what happens with raised-cosine QPSK: NRZ's near-constant envelope
    # concentrates a strong, prominent line at 4*fc, but continuous
    # envelope variation from pulse shaping spreads 4th-power energy into
    # sidebands and can leave no single dominant peak inside the window.
    band_mag = mag[band]
    peak_mag = mag[idx]
    typical_mag = np.median(band_mag)
    prominence_db = 20.0 * np.log10((peak_mag + 1e-300) / (typical_mag + 1e-300))
    if prominence_db < 6.0:
        print("  WARNING: carrier search found no clearly dominant peak "
              f"(only {prominence_db:.1f} dB above the typical level across "
              f"the +/-{search_frac*100:.0f}% window) -- this estimate is "
              "likely unreliable, not a genuine spectral line. Consider "
              "--fc-fixed with a known-good value from a companion "
              "capture on the same hardware.")

    # Quadratic interpolation on the log magnitude sharpens the peak
    if 0 < idx < len(mag) - 1:
        a, b, c = (np.log(mag[idx - 1] + 1e-300),
                   np.log(mag[idx] + 1e-300),
                   np.log(mag[idx + 1] + 1e-300))
        denom = a - 2 * b + c
        delta = 0.5 * (a - c) / denom if abs(denom) > 1e-30 else 0.0
        delta = float(np.clip(delta, -0.5, 0.5))
    else:
        delta = 0.0

    df = freqs[1] - freqs[0]
    f4 = freqs[idx] + delta * df

    return f4 / 4.0


# ===========================================================================
#  DEMODULATION
# ===========================================================================

def passband_to_baseband(x, fs, fc, rs, lpf_cutoff_hz, lpf_order=4):
    """
    Mix a real passband capture down to complex baseband.

    Multiplying by 2 * exp(-j*2*pi*fc*t) produces the wanted signal at DC
    and an unwanted image at -2*fc.  The low-pass filter removes the image.
    With fc = 1 kHz and Rs = 200 Bd the image sits at 2 kHz, ten times the
    symbol rate, so a modest filter is enough.

    filtfilt is used so the filter contributes no group delay, which keeps
    the symbol timing estimate honest.
    """
    n = len(x)
    t = np.arange(n) / fs

    y = x * (2.0 * np.exp(-1j * 2.0 * np.pi * fc * t))

    nyq = fs / 2.0
    wn = lpf_cutoff_hz / nyq
    if not (0.0 < wn < 1.0):
        raise ValueError(
            f"LPF cutoff {lpf_cutoff_hz:.1f} Hz is invalid for fs = {fs:.0f} Sa/s"
        )

    b, a = signal.butter(lpf_order, wn, btype="low")
    y = signal.filtfilt(b, a, y.real) + 1j * signal.filtfilt(b, a, y.imag)

    return y


def matched_filter(y, sps, shape):
    """
    Apply the receive filter matched to the transmit pulse.

    For rectangular NRZ the matched filter is an integrate-and-dump, which
    is a boxcar of one symbol length.  Convolving in 'same' mode keeps the
    output centred, so no delay compensation is needed.

    For raised-cosine shaping the transmit pulse already satisfies the
    Nyquist criterion on its own, so no further shaping is applied.  (Split
    the filtering between transmitter and receiver with root raised cosine
    if a true matched-filter pair is wanted.)
    """
    if shape == "nrz":
        h = np.ones(int(round(sps))) / float(int(round(sps)))
        return (np.convolve(y.real, h, mode="same")
                + 1j * np.convolve(y.imag, h, mode="same"))
    return y


# ===========================================================================
#  SYMBOL TIMING
# ===========================================================================

def oerder_meyr_offset(y, sps):
    """
    Oerder-Meyr square-law timing estimate.

    The squared magnitude of the filtered signal contains a spectral line
    at the symbol rate, produced by the amplitude dips at symbol
    transitions.  The phase of that line gives the timing offset directly,
    with no feedback loop.

    Returns a fractional offset in the range [0, 1) symbols.
    """
    n = len(y)
    k = np.arange(n)
    mag2 = np.abs(y) ** 2
    X = np.sum(mag2 * np.exp(-1j * 2.0 * np.pi * k / sps))
    eps = -np.angle(X) / (2.0 * np.pi)
    return float(eps % 1.0)


def extract_symbols(y, fs, rs, tau_frac, guard_sym=0.5):
    """
    Sample the filtered baseband signal at the symbol instants.

    Cubic spline interpolation is used because the symbol instants almost
    never land exactly on a sample: the scope clock and the STM32 clock are
    independent, so sps is not an integer in general.

    A guard of half a symbol is dropped at each end to avoid the
    matched-filter transient.
    """
    n = len(y)
    idx = np.arange(n, dtype=float)

    cs_r = CubicSpline(idx, y.real)
    cs_i = CubicSpline(idx, y.imag)

    sps = fs / rs
    start = guard_sym * sps + tau_frac * sps
    stop = n - 1 - guard_sym * sps

    if stop <= start:
        return np.array([], dtype=complex), np.array([])

    k_max = int(np.floor((stop - start) / sps))
    pos = start + np.arange(k_max + 1) * sps
    pos = pos[pos <= stop]

    z = cs_r(pos) + 1j * cs_i(pos)
    return z, pos / fs


def refine_timing(y, fs, rs, tau_init, n_grid=201, span=0.5):
    """
    Refine the timing offset by minimising EVM directly.

    Oerder-Meyr gives a good starting point but assumes a fairly long
    record.  With only eight to sixteen symbols available it can be a few
    percent of a symbol out, which matters at these EVM levels.  A direct
    search over the offset is cheap here and needs no assumptions about
    record length.
    """
    offsets = tau_init + np.linspace(-span, span, n_grid)
    best_tau, best_evm = tau_init, np.inf

    for tau in offsets:
        z, _ = extract_symbols(y, fs, rs, tau % 1.0)
        if len(z) < 4:
            continue
        zc, _, _ = align_phase(z, np.arange(len(z)) / rs)
        evm, _ = evm_rms(zc)
        if evm < best_evm:
            best_evm, best_tau = evm, tau % 1.0

    return best_tau, best_evm


# ===========================================================================
#  FREQUENCY, PHASE AND AMBIGUITY
# ===========================================================================

def align_phase(z, t_sym):
    """
    Remove residual frequency offset and static phase rotation.

    The fourth power trick applies to the recovered symbols too.  Raising
    them to the fourth power strips the modulation, leaving a slowly
    rotating phasor whose rotation rate is four times the residual carrier
    frequency error.  A straight line fitted to the unwrapped phase gives
    both the frequency error (slope) and the static rotation (intercept).

    Returns (corrected symbols, residual frequency in Hz, phase in radians).
    """
    if len(z) < 3:
        return z, 0.0, 0.0

    z4 = z ** 4
    phase = np.unwrap(np.angle(z4))

    slope, intercept = np.polyfit(t_sym, phase, 1)

    df = slope / (4.0 * 2.0 * np.pi)

    # Ideal points sit at 45 deg + m*90 deg, so their fourth power sits at
    # 180 deg.  Subtracting pi leaves the true static rotation.
    phi = (intercept - np.pi) / 4.0

    zc = z * np.exp(-1j * (2.0 * np.pi * df * t_sym + phi))

    return zc, float(df), float(phi)


def normalise_power(z):
    """Scale so the mean symbol power is 1, matching the ideal points."""
    p = np.mean(np.abs(z) ** 2)
    return z / np.sqrt(p) if p > 0 else z


def decide(z):
    """Hard-decide each symbol to the nearest ideal point."""
    d = np.zeros(len(z), dtype=int)
    d[z.real < 0] += 1     # b0 = 1
    d[z.imag < 0] += 2     # b1 = 1
    return d


def resolve_ambiguity(z, tx_pattern):
    """
    Resolve the 4-fold phase ambiguity and the unknown pattern start.

    A QPSK receiver cannot tell which of the four constellation points was
    meant to be 45 degrees, because all four look identical after carrier
    recovery.  In a real link a differential encoder or a known preamble
    settles this.  Here the transmitted pattern is known, so every
    combination of the four rotations and every cyclic start position is
    tried, and the one with fewest bit errors wins.

    Returns (rotated symbols, reference symbol sequence, bit errors, rotation).
    """
    best = None

    for m in range(4):
        zr = z * np.exp(1j * m * np.pi / 2.0)
        dec = decide(zr)

        for shift in range(len(tx_pattern)):
            ref = np.array([tx_pattern[(shift + k) % len(tx_pattern)]
                            for k in range(len(dec))], dtype=int)

            # Compare bit by bit, not symbol by symbol, so Gray coding
            # is credited properly.
            errs = int(np.sum(((dec & 1) != (ref & 1)).astype(int)
                              + (((dec >> 1) & 1) != ((ref >> 1) & 1)).astype(int)))

            if best is None or errs < best[2]:
                best = (zr, ref, errs, m)

    return best


# ===========================================================================
#  METRICS
# ===========================================================================

def evm_rms(z, ref_sym=None):
    """
    Error vector magnitude, in percent.

    EVM is the RMS distance between each received symbol and its ideal
    position, divided by the RMS ideal amplitude.  It is the standard
    single-number quality figure for a digitally modulated transmitter,
    and it folds phase error, amplitude error, noise and distortion into
    one value.
    """
    zn = normalise_power(z)
    ideal = IDEAL[decide(zn)] if ref_sym is None else IDEAL[ref_sym]
    err = zn - ideal
    num = np.mean(np.abs(err) ** 2)
    den = np.mean(np.abs(ideal) ** 2)
    return 100.0 * np.sqrt(num / den), err


def iq_imbalance(z):
    """
    Gain imbalance in dB and quadrature error in degrees.

    These matter most in baseband mode, where I and Q take separate paths
    through separate DAC channels and separate filter sections.  A gain
    mismatch stretches the constellation into a rectangle; a quadrature
    error shears it into a diamond.
    """
    i, q = z.real, z.imag
    si, sq = np.std(i), np.std(q)

    gain_db = 20.0 * np.log10(si / sq) if sq > 0 else np.nan

    if si > 0 and sq > 0:
        rho = np.mean(i * q) / (si * sq)
        quad_deg = np.degrees(np.arcsin(np.clip(rho, -1.0, 1.0)))
    else:
        quad_deg = np.nan

    return float(gain_db), float(quad_deg)


def phase_error_stats(z, ref_sym):
    """RMS and peak phase error against the ideal constellation, in degrees."""
    zn = normalise_power(z)
    ideal = IDEAL[ref_sym]
    err = np.angle(zn * np.conj(ideal))
    return float(np.degrees(np.sqrt(np.mean(err ** 2)))), \
           float(np.degrees(np.max(np.abs(err))))


# ===========================================================================
#  PLOTS
# ===========================================================================

def plot_waveform(t, x, prefix, label="Point B output"):
    fig, ax = plt.subplots(figsize=(7.0, 3.2))
    ax.plot(t * 1e3, x, lw=0.8, color="#1f4e79")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Voltage (V)")
    ax.set_title(f"QPSK capture, {label}")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{prefix}_waveform.png", dpi=300)
    plt.close(fig)


def plot_spectrum(x, fs, prefix, fc=None, rs=None):
    """
    Single-record spectrum with a Blackman-Harris window.

    Blackman-Harris is used throughout this project because the STM32 and
    the scope run from independent clocks, so no capture is ever coherent.
    Leakage, not resolution, is the limiting factor, and the -92 dB
    sidelobes of Blackman-Harris keep sidebands 50 to 60 dB down visible.
    Welch averaging is not available: a 2000-point record cannot be split
    into useful segments.
    """
    n = len(x)
    w = signal.windows.blackmanharris(n)
    cg = np.sum(w) / n

    X = np.fft.rfft((x - np.mean(x)) * w) / (n * cg)
    freqs = np.fft.rfftfreq(n, d=1.0 / fs)
    mag = 20.0 * np.log10(np.abs(X) + 1e-15)
    mag -= np.max(mag)

    fig, ax = plt.subplots(figsize=(7.0, 3.6))
    ax.plot(freqs, mag, lw=0.8, color="#1f4e79")

    if fc is not None:
        ax.axvline(fc, color="#c00000", ls="--", lw=0.8,
                   label=f"fc = {fc:.2f} Hz")
        if rs is not None:
            for s in (-1, 1):
                ax.axvline(fc + s * rs, color="#7f7f7f", ls=":", lw=0.8)
            ax.plot([], [], color="#7f7f7f", ls=":", lw=0.8,
                    label=f"fc +/- Rs ({rs:.2f} Bd)")
        ax.legend(fontsize=8)

    ax.set_xlim(0, min(5.0 * (fc or 1000.0), fs / 2.0))
    ax.set_ylim(-90, 5)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Magnitude (dB, normalised)")
    ax.set_title("QPSK spectrum (Blackman-Harris window)")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{prefix}_spectrum.png", dpi=300)
    plt.close(fig)


def plot_iq_tracks(y, fs, t_sym, z, prefix):
    """
    z is used only to fetch symbol timing context from the caller's
    workflow; the markers themselves are re-sampled directly from y at
    t_sym, not taken from the caller's (possibly frequency/phase
    corrected and power-normalised) symbol array. Plotting corrected,
    rescaled markers on top of the raw curve makes them visibly sit off
    the curve -- consistently, by exactly the correction and
    normalisation factors applied -- which looks like a fault even
    though both arrays are individually correct. Re-sampling y directly
    guarantees the markers sit exactly on the curve shown, since they
    come from the same data by construction.
    """
    idx = np.arange(len(y), dtype=float)
    cs_r = CubicSpline(idx, y.real)
    cs_i = CubicSpline(idx, y.imag)
    pos = t_sym * fs
    z_on_curve = cs_r(pos) + 1j * cs_i(pos)

    t = np.arange(len(y)) / fs
    fig, ax = plt.subplots(figsize=(7.0, 3.6))
    ax.plot(t * 1e3, y.real, lw=1.0, color="#1f4e79", label="I")
    ax.plot(t * 1e3, y.imag, lw=1.0, color="#c00000", label="Q")
    ax.plot(t_sym * 1e3, z_on_curve.real, "o", ms=4, color="#1f4e79", mfc="white")
    ax.plot(t_sym * 1e3, z_on_curve.imag, "o", ms=4, color="#c00000", mfc="white")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Amplitude (normalised)")
    ax.set_title("Recovered baseband I and Q, circles mark symbol instants")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{prefix}_iq_tracks.png", dpi=300)
    plt.close(fig)


def plot_constellation(z, ref_sym, evm, prefix):
    zn = normalise_power(z)
    ideal = IDEAL[ref_sym]

    fig, ax = plt.subplots(figsize=(4.8, 4.8))

    # Error vectors
    for a, b in zip(zn, ideal):
        ax.plot([b.real, a.real], [b.imag, a.imag],
                color="#c00000", lw=0.7, alpha=0.7, zorder=1)

    ax.plot(IDEAL.real, IDEAL.imag, "+", ms=14, mew=1.8,
            color="#444444", zorder=2, label="Ideal")
    ax.plot(zn.real, zn.imag, "o", ms=6, color="#1f4e79",
            mfc="#1f4e79", alpha=0.85, zorder=3, label="Measured")

    circ = plt.Circle((0, 0), 1.0, fill=False, color="#999999",
                      ls=":", lw=0.8)
    ax.add_patch(circ)

    ax.axhline(0, color="#bbbbbb", lw=0.6)
    ax.axvline(0, color="#bbbbbb", lw=0.6)
    ax.set_xlim(-1.6, 1.6)
    ax.set_ylim(-1.6, 1.6)
    ax.set_aspect("equal")
    ax.set_xlabel("In-phase (I)")
    ax.set_ylabel("Quadrature (Q)")
    ax.set_title(f"QPSK constellation, EVM = {evm:.2f} %")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(f"{prefix}_constellation.png", dpi=300)
    plt.close(fig)


def plot_eye(y, fs, rs, tau, prefix):
    """Two-symbol eye diagram of the in-phase track."""
    sps = fs / rs
    n = len(y)
    span = 2.0
    seg = int(round(span * sps))
    if seg < 8 or n < 3 * seg:
        return

    idx = np.arange(n, dtype=float)
    cs = CubicSpline(idx, y.real)

    fig, ax = plt.subplots(figsize=(5.6, 3.4))
    start = tau * sps
    k = 0
    while start + (k + span) * sps < n - 1:
        pos = start + k * sps + np.linspace(0, span * sps, seg)
        pos = pos[pos <= n - 1]
        ax.plot((pos - pos[0]) / sps, cs(pos), color="#1f4e79",
                lw=0.9, alpha=0.6)
        k += 1

    ax.set_xlabel("Time (symbol periods)")
    ax.set_ylabel("I amplitude")
    ax.set_title("Eye diagram, in-phase track")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{prefix}_eye.png", dpi=300)
    plt.close(fig)


# ===========================================================================
#  MAIN
# ===========================================================================

def process_capture(t, channels, args):
    """Turn one capture into a set of recovered symbols."""
    fs = 1.0 / np.mean(np.diff(t))

    if args.mode == "passband":
        x = channels[args.ch - 1]
        x = x - np.mean(x)

        fc = args.fc if args.fc_fixed else estimate_carrier(x, fs, args.fc)

        # The carrier and the symbol clock both come from TIM6, so they
        # share one clock error exactly.  Estimating fc therefore also
        # calibrates Rs, with no separate measurement needed.
        clock_ratio = fc / args.fc
        rs = args.rs * clock_ratio

        cutoff = args.lpf if args.lpf else 3.0 * rs
        y = passband_to_baseband(x, fs, fc, rs, cutoff)
        y = matched_filter(y, fs / rs, args.shape)

    else:
        if len(channels) < 2:
            raise ValueError(
                "Baseband mode needs two channels. Supply a two-channel "
                "CSV, or pass --csv-q with the Q-channel file."
            )
        i_ch = channels[0] - np.mean(channels[0])
        q_ch = channels[1] - np.mean(channels[1])

        fc = 0.0
        clock_ratio = 1.0
        rs = args.rs

        y = i_ch + 1j * q_ch
        if args.lpf:
            b, a = signal.butter(4, args.lpf / (fs / 2.0), btype="low")
            y = signal.filtfilt(b, a, y.real) + 1j * signal.filtfilt(b, a, y.imag)
        y = matched_filter(y, fs / rs, args.shape)

    tau0 = oerder_meyr_offset(y, fs / rs)
    tau, _ = refine_timing(y, fs, rs, tau0)

    z_raw, t_sym = extract_symbols(y, fs, rs, tau)
    if len(z_raw) < 4:
        raise ValueError(
            f"Only {len(z_raw)} symbols recovered. Use a slower timebase "
            "so the record spans more symbols."
        )

    z, df, phi = align_phase(z_raw, t_sym)
    z = normalise_power(z)

    # Resolve rotation and cyclic pattern shift for THIS capture only.
    #
    # Each single-shot trigger grabs the waveform at an arbitrary point in
    # the continuously-repeating symbol pattern, and each capture's own
    # 4th-power phase fit lands on one of four equally valid 90-degree
    # branches independently of any other capture. Neither the pattern
    # start position nor the rotation branch is shared between separate
    # trigger events. Resolving both per capture -- rather than once over
    # a multi-capture pool -- is what makes pooling several short captures
    # valid at all.
    zr, ref, bit_errs, rot = resolve_ambiguity(z, TX_PATTERN)
    zr = normalise_power(zr)

    return {
        "fs": fs, "fc": fc, "rs": rs, "clock_ratio": clock_ratio,
        "tau": tau, "tau_om": tau0, "df": df, "phi": phi,
        "y": y, "z": z, "t_sym": t_sym, "x": channels[args.ch - 1], "t": t,
        "z_corrected": zr, "ref": ref, "bit_errs": bit_errs,
        "n_bits": 2 * len(zr), "rotation_deg": rot * 90,
    }


def plot_snr_sweep(snr_db, evm_meas, evm_std, evm_theory, ber_meas, ber_theory, prefix):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.0, 4.2))

    ax1.errorbar(snr_db, evm_meas, yerr=evm_std, fmt="o-", ms=4, lw=1.2,
                 color="#1f4e79", label="Measured (mean +/- std)")
    ax1.plot(snr_db, evm_theory, "--", lw=1.2, color="#c00000",
              label="Theory, EVM ~ 1/sqrt(SNR)")
    ax1.set_xlabel("SNR (dB)")
    ax1.set_ylabel("EVM (%)")
    ax1.set_yscale("log")
    ax1.set_title("EVM against SNR")
    ax1.grid(True, which="both", alpha=0.3)
    ax1.legend(fontsize=8)

    ax2.semilogy(snr_db, np.clip(ber_meas, 1e-4, None), "o-", ms=4, lw=1.2,
                 color="#1f4e79", label="Measured")
    ax2.semilogy(snr_db, np.clip(ber_theory, 1e-12, None), "--", lw=1.2,
                 color="#c00000", label="Theory, Q(sqrt(2*SNR))")
    ax2.set_xlabel("SNR (dB)")
    ax2.set_ylabel("BER")
    ax2.set_title("QPSK BER against SNR (approximate, SNR taken as Eb/N0)")
    ax2.grid(True, which="both", alpha=0.3)
    ax2.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(f"{prefix}_snr_sweep.png", dpi=300)
    plt.close(fig)


def run_snr_sweep(x_clean, fs, fc, rs, args):
    """
    EVM-vs-SNR sweep using synthetic noise added to one clean, measured
    QPSK capture. This explores parameter space the bench data alone
    cannot reach -- a software experiment on real transmitter data, not a
    new physical measurement. fc and rs are held fixed at the values
    already established from the clean capture, since carrier estimation
    on a heavily noised signal is itself unreliable at low SNR.

    SNR is defined as signal power to added-noise power over the full
    passband record, used as an approximation to Eb/N0 for the
    theoretical BER curve. The two are not identical in general -- the
    exact relationship depends on pulse shape, samples per symbol, and
    receiver filtering -- so this is stated as an approximation, not an
    exact equivalence. QPSK's theoretical BER uses the same Q-function as
    BPSK per bit, since Gray coding makes each bit's error probability
    independent of the other bit in the same symbol.
    """
    from scipy.special import erfc

    rng = np.random.default_rng(args.snr_seed)
    p_signal = np.mean(x_clean ** 2)

    snr_range = np.arange(args.snr_min, args.snr_max + 1e-9, args.snr_step)
    n_real = args.snr_realizations

    evm_mean, evm_std, ber_mean = [], [], []

    print()
    print("-" * 68)
    print("EVM-VS-SNR SWEEP")
    print("-" * 68)
    print(f"Clean capture power    : {p_signal:.6g} (arb. units)")
    print(f"Realizations per point : {n_real}")
    print(f"SNR range              : {args.snr_min:.1f} to {args.snr_max:.1f} dB, "
          f"step {args.snr_step:.1f} dB")
    print()
    print(f"{'SNR (dB)':>10}  {'EVM mean (%)':>13}  {'EVM std (%)':>12}  "
          f"{'BER':>10}")

    for snr_db in snr_range:
        noise_power = p_signal / (10.0 ** (snr_db / 10.0))
        evms, errs_total, bits_total = [], 0, 0

        for _ in range(n_real):
            noise = rng.normal(0.0, np.sqrt(noise_power), size=x_clean.shape)
            x_noisy = x_clean + noise

            cutoff = args.lpf if args.lpf else 3.0 * rs
            y = passband_to_baseband(x_noisy, fs, fc, rs, cutoff)
            y = matched_filter(y, fs / rs, args.shape)

            tau0 = oerder_meyr_offset(y, fs / rs)
            tau, _ = refine_timing(y, fs, rs, tau0)
            z_raw, t_sym = extract_symbols(y, fs, rs, tau)
            if len(z_raw) < 4:
                continue

            z, _, _ = align_phase(z_raw, t_sym)
            z = normalise_power(z)
            zr, ref, bit_errs, _ = resolve_ambiguity(z, TX_PATTERN)
            zr = normalise_power(zr)

            evm, _ = evm_rms(zr, ref)
            evms.append(evm)
            errs_total += bit_errs
            bits_total += 2 * len(zr)

        if not evms:
            continue

        evm_mean.append(np.mean(evms))
        evm_std.append(np.std(evms))
        ber_mean.append(errs_total / bits_total if bits_total else np.nan)

        print(f"{snr_db:>10.1f}  {evm_mean[-1]:>13.3f}  {evm_std[-1]:>12.3f}  "
              f"{ber_mean[-1]:>10.3e}")

    snr_db_arr = snr_range[:len(evm_mean)]
    evm_mean = np.array(evm_mean)
    evm_std = np.array(evm_std)
    ber_mean = np.array(ber_mean)

    snr_lin = 10.0 ** (snr_db_arr / 10.0)
    evm_theory = 100.0 / np.sqrt(snr_lin)
    ber_theory = 0.5 * erfc(np.sqrt(snr_lin))

    plot_snr_sweep(snr_db_arr, evm_mean, evm_std, evm_theory,
                   ber_mean, ber_theory, args.prefix)

    with open(f"{args.prefix}_snr_sweep.txt", "w", encoding="utf-8") as f:
        f.write("EVM-vs-SNR sweep, synthetic noise on a measured QPSK capture\n")
        f.write("=" * 62 + "\n")
        f.write(f"Realizations per point: {n_real}\n\n")
        f.write(f"{'SNR(dB)':>10} {'EVM_meas(%)':>13} {'EVM_std(%)':>12} "
                f"{'EVM_theory(%)':>14} {'BER_meas':>12} {'BER_theory':>12}\n")
        for i in range(len(snr_db_arr)):
            f.write(f"{snr_db_arr[i]:>10.1f} {evm_mean[i]:>13.3f} "
                    f"{evm_std[i]:>12.3f} {evm_theory[i]:>14.3f} "
                    f"{ber_mean[i]:>12.3e} {ber_theory[i]:>12.3e}\n")

    print()
    print(f"  {args.prefix}_snr_sweep.png")
    print(f"  {args.prefix}_snr_sweep.txt")


def main():
    p = argparse.ArgumentParser(
        description="QPSK constellation recovery and EVM analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument("--csv", nargs="+", required=True,
                   help="One or more scope CSV files. Several files are "
                        "pooled to improve the EVM statistics.")
    p.add_argument("--csv-q", nargs="+", default=None,
                   help="Q-channel files, if I and Q were saved separately.")
    p.add_argument("--mode", choices=["passband", "baseband"],
                   default="passband",
                   help="passband = firmware MODE 0, one channel on PA4. "
                        "baseband = firmware MODE 1, I on PA4 and Q on PA5.")
    p.add_argument("--fc", type=float, default=FC_NOMINAL_HZ,
                   help="Nominal carrier frequency in Hz (default 1000).")
    p.add_argument("--fc-fixed", action="store_true",
                   help="Trust --fc exactly instead of estimating it.")
    p.add_argument("--rs", type=float, default=RS_NOMINAL_HZ,
                   help="Nominal symbol rate in Bd (default 200).")
    p.add_argument("--shape", choices=["nrz", "rc"], default="nrz",
                   help="Transmit pulse shape, must match PULSE_SHAPE "
                        "in the firmware.")
    p.add_argument("--lpf", type=float, default=None,
                   help="Low-pass cutoff in Hz (default 3 x Rs).")
    p.add_argument("--ch", type=int, default=1,
                   help="Which channel holds the signal in passband mode.")
    p.add_argument("--prefix", default="qpsk",
                   help="Output filename prefix.")
    p.add_argument("--no-plots", action="store_true")

    p.add_argument("--snr-sweep", action="store_true",
                   help="Run an EVM-vs-SNR sweep using the first --csv "
                        "file as the noise-free basis.")
    p.add_argument("--snr-min", type=float, default=-5.0)
    p.add_argument("--snr-max", type=float, default=20.0)
    p.add_argument("--snr-step", type=float, default=1.0)
    p.add_argument("--snr-realizations", type=int, default=30,
                   help="Independent noise draws averaged per SNR point.")
    p.add_argument("--snr-seed", type=int, default=0)

    args = p.parse_args()

    print("=" * 68)
    print("QPSK CONSTELLATION RECOVERY AND EVM ANALYSIS")
    print("=" * 68)
    print(f"Mode          : {args.mode}")
    print(f"Pulse shape   : {args.shape}")
    print(f"Nominal fc    : {args.fc:.1f} Hz")
    print(f"Nominal Rs    : {args.rs:.1f} Bd")
    print()
    print("Loading captures")

    results = []
    for k, path in enumerate(args.csv):
        t, chans = load_tek_csv(path)

        if args.csv_q is not None and k < len(args.csv_q):
            _, chans_q = load_tek_csv(args.csv_q[k])
            chans = [chans[0], chans_q[0]]

        try:
            r = process_capture(t, chans, args)
            r["file"] = os.path.basename(path)
            results.append(r)
        except Exception as exc:
            print(f"  !! {os.path.basename(path)} skipped: {exc}")

    if not results:
        print("\nNo capture could be processed.")
        sys.exit(1)

    print()
    print("-" * 68)
    print("PER-CAPTURE RESOLUTION")
    print("-" * 68)
    print("Rotation and cyclic pattern shift are resolved independently for")
    print("each capture, since separate trigger events share neither an")
    print("absolute carrier phase reference nor a common position in the")
    print("repeating symbol pattern.")
    print()
    for r in results:
        n_sym = len(r["z_corrected"])
        print(f"  {r['file']}: rotation={r['rotation_deg']:>3.0f} deg, "
              f"symbols={n_sym}, bit errors={r['bit_errs']} of {r['n_bits']}")

    # Pool the per-capture CORRECTED symbols and their matched reference
    # sequences. Each capture has already been independently rotated and
    # aligned to the known pattern above, so this concatenation combines
    # captures that are all expressed in the same absolute constellation
    # frame -- unlike pooling the raw symbols and resolving once.
    z_all = np.concatenate([r["z_corrected"] for r in results])
    ref_all = np.concatenate([r["ref"] for r in results])
    z_all = normalise_power(z_all)

    evm, err_vec = evm_rms(z_all, ref_all)
    evm_peak = 100.0 * np.max(np.abs(err_vec)) / np.sqrt(np.mean(np.abs(IDEAL[ref_all]) ** 2))
    mer_db = -20.0 * np.log10(evm / 100.0)

    gain_db, quad_deg = iq_imbalance(z_all)
    ph_rms, ph_pk = phase_error_stats(z_all, ref_all)

    bit_errs = int(sum(r["bit_errs"] for r in results))
    n_bits = int(sum(r["n_bits"] for r in results))

    r0 = results[0]
    zr = z_all

    print()
    print("-" * 68)
    print("RECOVERY")
    print("-" * 68)
    if args.mode == "passband":
        ppm = (r0["clock_ratio"] - 1.0) * 1e6
        print(f"Carrier estimate      : {r0['fc']:.4f} Hz "
              f"({ppm:+.0f} ppm vs nominal)")
        print(f"Symbol rate (scaled)  : {r0['rs']:.4f} Bd")
    print(f"Scope sample rate     : {r0['fs'] / 1e3:.3f} kSa/s")
    print(f"Samples per symbol    : {r0['fs'] / r0['rs']:.2f}")
    print(f"Timing offset (O-M)   : {r0['tau_om']:.4f} symbol")
    print(f"Timing offset (final) : {r0['tau']:.4f} symbol")
    print(f"Residual freq offset  : {r0['df']:+.4f} Hz")
    print(f"Static phase rotation : {np.degrees(r0['phi']):+.2f} deg")
    print(f"Ambiguity resolution   : per capture, see table above")

    print()
    print("-" * 68)
    print("RESULTS")
    print("-" * 68)
    print(f"Captures pooled       : {len(results)}")
    print(f"Symbols analysed      : {len(zr)}  ({n_bits} bits)")
    print(f"EVM (RMS)             : {evm:.3f} %")
    print(f"EVM (peak)            : {evm_peak:.3f} %")
    print(f"MER                   : {mer_db:.2f} dB")
    print(f"Phase error, RMS      : {ph_rms:.3f} deg")
    print(f"Phase error, peak     : {ph_pk:.3f} deg")
    print(f"I/Q gain imbalance    : {gain_db:+.3f} dB")
    print(f"Quadrature error      : {quad_deg:+.3f} deg")
    print(f"Bit errors            : {bit_errs} of {n_bits}")
    print(f"BER                   : {bit_errs / n_bits:.3e}")

    print()
    if evm < 1.0:
        print("NOTE: EVM is below the 8-bit instrument floor of roughly")
        print("0.6 percent. Report this as an upper bound on transmitter")
        print("EVM, not as a measured transmitter figure.")
    else:
        print("EVM is above the estimated 8-bit instrument floor of roughly")
        print("0.6 percent, so the result reflects the transmitter and its")
        print("analogue path rather than the digitiser alone.")

    if not args.no_plots:
        print()
        print("Writing figures")
        plot_waveform(r0["t"], r0["x"], args.prefix)
        plot_spectrum(r0["x"] - np.mean(r0["x"]), r0["fs"], args.prefix,
                      fc=(r0["fc"] if args.mode == "passband" else None),
                      rs=r0["rs"])
        plot_iq_tracks(r0["y"], r0["fs"], r0["t_sym"], r0["z"], args.prefix)
        plot_constellation(zr, ref_all, evm, args.prefix)
        plot_eye(r0["y"], r0["fs"], r0["rs"], r0["tau"], args.prefix)
        for suffix in ("waveform", "spectrum", "iq_tracks",
                       "constellation", "eye"):
            f = f"{args.prefix}_{suffix}.png"
            if os.path.isfile(f):
                print(f"  {f}")

    with open(f"{args.prefix}_results.txt", "w", encoding="utf-8") as f:
        f.write("QPSK constellation recovery and EVM analysis\n")
        f.write("=" * 50 + "\n")
        f.write(f"Mode                  : {args.mode}\n")
        f.write(f"Pulse shape           : {args.shape}\n")
        f.write(f"Files                 : {', '.join(args.csv)}\n")
        f.write(f"Carrier estimate      : {r0['fc']:.4f} Hz\n")
        f.write(f"Symbol rate           : {r0['rs']:.4f} Bd\n")
        f.write(f"Scope sample rate     : {r0['fs']:.1f} Sa/s\n")
        f.write(f"Symbols analysed      : {len(zr)}\n")
        f.write(f"EVM (RMS)             : {evm:.3f} %\n")
        f.write(f"EVM (peak)            : {evm_peak:.3f} %\n")
        f.write(f"MER                   : {mer_db:.2f} dB\n")
        f.write(f"Phase error, RMS      : {ph_rms:.3f} deg\n")
        f.write(f"I/Q gain imbalance    : {gain_db:+.3f} dB\n")
        f.write(f"Quadrature error      : {quad_deg:+.3f} deg\n")
        f.write(f"Bit errors            : {bit_errs} of {n_bits}\n")
    print(f"  {args.prefix}_results.txt")

    if args.snr_sweep:
        run_snr_sweep(r0["x"], r0["fs"], r0["fc"], r0["rs"], args)


if __name__ == "__main__":
    main()
