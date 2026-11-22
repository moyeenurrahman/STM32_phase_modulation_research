#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bpsk_analysis.py
================

BPSK software-defined receiver and EVM analysis for STM32 baseband captures.

Project : Implementation and Performance Analysis of a Baseband Phase
          Modulation System with STM32 for Low-Power Wireless Communication
Board   : NUCLEO-F303RE, DAC1 + TIM6 TRGO + circular DMA
Scope   : Tektronix TBS1000C, 2000-point record, 8-bit vertical

This is the receive-side companion to the transmit-side QPSK work
(qpsk_analysis.py). It demodulates a real, single-channel BPSK capture
entirely in software -- carrier recovery, matched filtering, symbol timing,
phase alignment, decision -- with no assumptions beyond the nominal carrier
and symbol rate. This is the "software-defined radio" half of the story
from the receive side: a generic sampled waveform in, a decoded bit stream
and quality metrics out, with no hardware specific to BPSK involved beyond
the ADC-like samples themselves.

Processing chain
-----------------
  1. Load and de-trend the capture
  2. Estimate the true carrier frequency by the squaring method
  3. Scale the symbol rate by the same clock error
  4. Mix down to complex baseband, low-pass filter
  5. Matched filter (integrate-and-dump for NRZ, pass-through for RC)
  6. Symbol timing recovery: Oerder-Meyr, refined by minimum-EVM search
  7. Residual frequency and phase correction (squaring loop)
  8. Resolve the 2-fold phase ambiguity against the known bit pattern
  9. Compute EVM, MER, phase error, and bit errors

Typical use
-----------
  # Single capture
  python bpsk_analysis.py --csv capture.csv --fc 5000 --rs 500 \\
      --pattern 0,1,0,1,0,1,0,1

  # Several captures pooled for better EVM statistics
  python bpsk_analysis.py --csv cap1.csv cap2.csv cap3.csv \\
      --fc 5000 --rs 500 --pattern 0,1,0,1,0,1,0,1

  # Raised-cosine shaped capture
  python bpsk_analysis.py --csv capture.csv --shape rc --fc 5000 --rs 500 \\
      --pattern 0,1,0,1,0,1,0,1

  # EVM-vs-SNR sweep, using one clean capture as the noise-free basis
  python bpsk_analysis.py --csv clean_capture.csv --fc 5000 --rs 500 \\
      --pattern 0,1,0,1,0,1,0,1 --snr-sweep

IMPORTANT -- confirm before running on real data
-------------------------------------------------
--fc, --rs and --pattern MUST match the actual firmware that produced the
capture, not the defaults below. The defaults are placeholders only. Running
this against the wrong pattern or symbol rate does not fail loudly -- it
silently produces a plausible-looking but meaningless EVM and BER, exactly
the failure mode this project has already hit twice while building the
QPSK companion script. When in doubt, open the firmware's symbol/bit array
and transcribe it into --pattern exactly.
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
#  DEFAULTS -- PLACEHOLDERS, CONFIRM AGAINST YOUR FIRMWARE BEFORE TRUSTING
# ===========================================================================

FC_NOMINAL_HZ = 5000.0    # confirmed from phase-accuracy BPSK main.c
RS_NOMINAL_HZ = 250.0     # confirmed from phase-accuracy BPSK main.c

# BPSK constellation: bit 0 -> +1 (0 deg), bit 1 -> -1 (180 deg).
# If your firmware uses the opposite mapping, every bit comes out inverted
# but EVM and BER are unaffected (both bits are still compared symmetrically
# by resolve_ambiguity's 2-way rotation search).
IDEAL = np.array([+1.0, -1.0])

# Confirmed transmitted bit pattern, from the phase-accuracy BPSK main.c
# (static const uint8_t bits[10] = {0,1,0,1,1,0,1,0,0,1};). Override with
# --pattern if a different firmware build produced the capture being
# analysed.
TX_PATTERN_DEFAULT = np.array([0, 1, 0, 1, 1, 0, 1, 0, 0, 1], dtype=int)


# ===========================================================================
#  CSV LOADING (shared logic with qpsk_analysis.py)
# ===========================================================================

def _to_float(s):
    try:
        return float(s.strip())
    except (ValueError, AttributeError):
        return np.nan


def load_tek_csv(path, verbose=True):
    """
    Load a Tektronix CSV export and return (t, [ch1, ch2, ...]).

    Format-agnostic: locates the longest run of rows whose right-hand
    fields all parse as numbers, then identifies the time axis as the
    column with a positive, near-constant step.

    Tolerance note: OpenChoice CSV exports write timestamps as text with
    about 9 significant figures, which alone produces step-to-step
    relative variation on the order of 1e-5. A tolerance tighter than
    that rejects every real capture; 1e-2 (1 percent) accepts normal
    export rounding while still rejecting a genuinely oscillating signal
    column, which fails the strict positivity check regardless.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(path)

    with open(path, "r", newline="", encoding="utf-8", errors="replace") as f:
        rows = [r for r in csv.reader(f)]

    if not rows:
        raise ValueError(f"{path}: file is empty")

    ncols = max(len(r) for r in rows)
    mat = np.full((len(rows), ncols), np.nan)
    for i, r in enumerate(rows):
        for j, cell in enumerate(r):
            mat[i, j] = _to_float(cell)

    numeric_mask = ~np.isnan(mat)
    signatures = [tuple(row) for row in numeric_mask]

    best_start, best_len = 0, 0
    i = 0
    while i < len(signatures):
        j = i
        while j + 1 < len(signatures) and signatures[j + 1] == signatures[i]:
            j += 1
        run = j - i + 1
        if sum(signatures[i]) >= 2 and run > best_len:
            best_start, best_len = i, run
        i = j + 1

    if best_len < 10:
        raise ValueError(f"{path}: could not find a waveform data block")

    block = mat[best_start:best_start + best_len, :]
    cols = np.where(~np.isnan(block[0]))[0]
    data = block[:, cols]

    time_col = None
    for k in range(data.shape[1]):
        col = data[:, k]
        d = np.diff(col)
        if np.all(d > 0) and np.std(d) < 1e-2 * abs(np.mean(d)) + 1e-15:
            time_col = k
            break

    if time_col is None:
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
#  CARRIER ESTIMATION -- SQUARING METHOD (BPSK, not 4th power)
# ===========================================================================

def estimate_carrier(x, fs, fc_nominal, search_frac=0.05):
    """
    Estimate the true carrier frequency by the squaring method.

    A BPSK signal carries no discrete line at fc, because its two phase
    states (0 and 180 degrees) average to zero. Squaring the signal
    collapses both states onto the same phase, because 2*0 deg and
    2*180 deg are both 0 mod 360, leaving a clean tone at 2*fc. This is
    the BPSK analogue of the 4th-power trick used for QPSK -- the power
    matches the number of phase states being collapsed.

    Warns if the located peak has poor prominence relative to the rest of
    the search band, which indicates no genuine spectral line was found
    (the estimate is then unreliable, not merely imprecise) rather than
    checking literal bin position -- a wide-main-lobe window can leak a
    local maximum several bins inside the window edge even when the true
    line lies outside it.
    """
    n = len(x)
    w = signal.windows.blackmanharris(n)
    x2 = (x * w) ** 2

    nfft = int(2 ** np.ceil(np.log2(n * 16)))
    X = np.fft.rfft(x2 - np.mean(x2), n=nfft)
    freqs = np.fft.rfftfreq(nfft, d=1.0 / fs)
    mag = np.abs(X)

    f_target = 2.0 * fc_nominal
    lo = f_target * (1.0 - search_frac)
    hi = f_target * (1.0 + search_frac)

    band = (freqs >= lo) & (freqs <= hi)
    if not np.any(band):
        raise ValueError(
            f"2*fc = {f_target:.1f} Hz falls outside the captured bandwidth "
            f"(fs = {fs:.0f} Sa/s). Use a faster timebase."
        )

    idx = np.argmax(np.where(band, mag, 0.0))

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
    f2 = freqs[idx] + delta * df

    return f2 / 2.0


# ===========================================================================
#  DEMODULATION
# ===========================================================================

def passband_to_baseband(x, fs, fc, rs, lpf_cutoff_hz, lpf_order=4):
    """Mix a real passband capture down to complex baseband and low-pass filter."""
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
    Integrate-and-dump for NRZ. For raised-cosine, the transmit pulse
    already satisfies the Nyquist zero-ISI criterion at the symbol
    centres, so no further receive-side shaping is applied.
    """
    if shape == "nrz":
        h = np.ones(int(round(sps))) / float(int(round(sps)))
        return (np.convolve(y.real, h, mode="same")
                + 1j * np.convolve(y.imag, h, mode="same"))
    return y


# ===========================================================================
#  SYMBOL TIMING (generic -- not specific to the number of phase states)
# ===========================================================================

def oerder_meyr_offset(y, sps):
    """Oerder-Meyr square-law timing estimate. Returns offset in [0, 1) symbols."""
    n = len(y)
    k = np.arange(n)
    mag2 = np.abs(y) ** 2
    X = np.sum(mag2 * np.exp(-1j * 2.0 * np.pi * k / sps))
    eps = -np.angle(X) / (2.0 * np.pi)
    return float(eps % 1.0)


def extract_symbols(y, fs, rs, tau_frac, guard_sym=0.5):
    """Sample the filtered baseband signal at the symbol instants (cubic spline)."""
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
    """Refine the timing offset by minimising EVM directly (short-record friendly).

    n_grid controls search resolution. The default (201) is used for the
    primary single-capture analysis, where timing precision matters most.
    A coarser grid is appropriate inside the SNR sweep, where the same
    search runs once per noise realization (potentially thousands of
    times) and noise-averaging already dominates the achievable precision
    far more than a finer timing grid would.
    """
    offsets = tau_init + np.linspace(-span, span, n_grid)
    best_tau, best_evm = tau_init, np.inf

    for tau in offsets:
        z, _ = extract_symbols(y, fs, rs, tau % 1.0)
        if len(z) < 3:
            continue
        zc, _, _ = align_phase(z, np.arange(len(z)) / rs)
        evm, _ = evm_rms(zc)
        if evm < best_evm:
            best_evm, best_tau = evm, tau % 1.0

    return best_tau, best_evm


# ===========================================================================
#  FREQUENCY, PHASE AND AMBIGUITY -- SQUARING LOOP (BPSK, not 4th power)
# ===========================================================================

def align_phase(z, t_sym):
    """
    Remove residual frequency offset and static phase rotation.

    Squaring the recovered symbols strips the modulation for BPSK: with
    ideal phases at 0 and 180 degrees, 2*0 = 0 and 2*180 = 360 = 0 mod
    360 -- both collapse to angle 0 (unlike QPSK's 4th-power collapse to
    180 degrees, since QPSK's phases are of the form (2k+1)*45 degrees).
    A line fitted to the unwrapped phase of z**2 gives both the residual
    frequency error (slope) and the static rotation (intercept).
    """
    if len(z) < 3:
        return z, 0.0, 0.0

    z2 = z ** 2
    phase = np.unwrap(np.angle(z2))

    slope, intercept = np.polyfit(t_sym, phase, 1)

    df = slope / (2.0 * 2.0 * np.pi)
    phi = intercept / 2.0

    zc = z * np.exp(-1j * (2.0 * np.pi * df * t_sym + phi))

    return zc, float(df), float(phi)


def normalise_power(z):
    p = np.mean(np.abs(z) ** 2)
    return z / np.sqrt(p) if p > 0 else z


def decide(z):
    """Hard-decide each symbol from the sign of the real part only."""
    return (z.real < 0).astype(int)


def resolve_ambiguity(z, tx_pattern):
    """
    Resolve the 2-fold phase ambiguity and the unknown pattern start.

    BPSK carrier recovery leaves a 2-fold (0/180 degree) ambiguity, and a
    single-shot capture starts at an arbitrary point in the repeating
    pattern. Both are resolved together by trying every rotation and
    every cyclic shift and keeping the combination with fewest bit
    errors against the known transmitted pattern.

    Also reports how many of the 16 (rotation, shift) combinations tie
    for the best error count. With very few symbols per capture -- as
    happens when the symbol period is long relative to the fixed scope
    record length -- a handful of decided bits can match multiple
    candidate alignments equally well by chance. A tie count above 1
    means this capture's resolution is not uniquely determined by its
    own data and should be treated cautiously, or excluded, rather than
    silently accepted.
    """
    best = None
    best_errs = None
    candidates = []

    for m in range(2):
        zr = z * np.exp(1j * m * np.pi)
        dec = decide(zr)

        for shift in range(len(tx_pattern)):
            ref = np.array([tx_pattern[(shift + k) % len(tx_pattern)]
                            for k in range(len(dec))], dtype=int)
            errs = int(np.sum(dec != ref))
            candidates.append(errs)

            if best is None or errs < best[2]:
                best = (zr, ref, errs, m)

    n_tied = sum(1 for e in candidates if e == best[2])
    if n_tied > 1 and len(z) <= 4:
        print(f"  WARNING: {n_tied} of {len(candidates)} rotation/shift "
              f"combinations tie for the best fit ({best[2]} errors) on "
              f"only {len(z)} symbols -- this capture's alignment is not "
              "uniquely determined and may be unreliable. Consider "
              "excluding it or pooling many more captures.")

    return (*best, n_tied)


# ===========================================================================
#  METRICS
# ===========================================================================

def evm_rms(z, ref_sym=None):
    """EVM in percent: RMS error-vector magnitude relative to the RMS ideal amplitude."""
    zn = normalise_power(z)
    ideal = IDEAL[decide(zn)] if ref_sym is None else IDEAL[ref_sym]
    err = zn - ideal
    num = np.mean(np.abs(err) ** 2)
    den = np.mean(np.abs(ideal) ** 2)
    return 100.0 * np.sqrt(num / den), err


def phase_error_stats(z, ref_sym):
    """RMS and peak phase error against the ideal constellation, in degrees."""
    zn = normalise_power(z)
    ideal = IDEAL[ref_sym]
    err = np.angle(zn * np.conj(ideal))
    return float(np.degrees(np.sqrt(np.mean(err ** 2)))), \
           float(np.degrees(np.max(np.abs(err))))


def residual_quadrature(z):
    """
    RMS of the imaginary part after normalisation, as a fraction of the
    total RMS amplitude -- diagnostic only. Ideal BPSK carries no
    information on Q; a large value here indicates leftover carrier
    phase error or an imperfect real-axis alignment, not a fault of BPSK
    itself.
    """
    zn = normalise_power(z)
    return float(100.0 * np.sqrt(np.mean(zn.imag ** 2)))


# ===========================================================================
#  PLOTS
# ===========================================================================

def plot_waveform(t, x, prefix):
    fig, ax = plt.subplots(figsize=(7.0, 3.2))
    ax.plot(t * 1e3, x, lw=0.8, color="#1f4e79")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Voltage (V)")
    ax.set_title("BPSK capture, Point B output")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{prefix}_waveform.png", dpi=300)
    plt.close(fig)


def plot_spectrum(x, fs, prefix, fc=None, rs=None):
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
        ax.axvline(fc, color="#c00000", ls="--", lw=0.8, label=f"fc = {fc:.2f} Hz")
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
    ax.set_title("BPSK spectrum (Blackman-Harris window)")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{prefix}_spectrum.png", dpi=300)
    plt.close(fig)


def plot_iq_tracks(y, fs, t_sym, z, prefix):
    """
    Markers are re-sampled directly from y at t_sym rather than taken
    from the caller's frequency/phase-corrected, power-normalised z, so
    they sit exactly on the curve shown. See qpsk_analysis.py's
    plot_iq_tracks for the full explanation of why this matters.
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

    for a, b in zip(zn, ideal):
        ax.plot([b, a.real], [0, a.imag], color="#c00000", lw=0.7,
                alpha=0.7, zorder=1)

    ax.plot(IDEAL, [0, 0], "+", ms=14, mew=1.8, color="#444444",
            zorder=2, label="Ideal")
    ax.plot(zn.real, zn.imag, "o", ms=6, color="#1f4e79", mfc="#1f4e79",
            alpha=0.85, zorder=3, label="Measured")

    circ = plt.Circle((0, 0), 1.0, fill=False, color="#999999", ls=":", lw=0.8)
    ax.add_patch(circ)

    ax.axhline(0, color="#bbbbbb", lw=0.6)
    ax.axvline(0, color="#bbbbbb", lw=0.6)
    ax.set_xlim(-1.6, 1.6)
    ax.set_ylim(-1.6, 1.6)
    ax.set_aspect("equal")
    ax.set_xlabel("In-phase (I)")
    ax.set_ylabel("Quadrature (Q)")
    ax.set_title(f"BPSK constellation, EVM = {evm:.2f} %")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(f"{prefix}_constellation.png", dpi=300)
    plt.close(fig)


def plot_eye(y, fs, rs, tau, prefix):
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
        ax.plot((pos - pos[0]) / sps, cs(pos), color="#1f4e79", lw=0.9, alpha=0.6)
        k += 1

    ax.set_xlabel("Time (symbol periods)")
    ax.set_ylabel("I amplitude")
    ax.set_title("Eye diagram, in-phase track")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{prefix}_eye.png", dpi=300)
    plt.close(fig)


def plot_snr_sweep(snr_db, evm_meas, evm_std, evm_theory_naive, evm_theory_mf,
                    ber_meas, ber_theory_naive, ber_theory_mf, prefix):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.0, 4.2))

    ax1.errorbar(snr_db, evm_meas, yerr=evm_std, fmt="o-", ms=4, lw=1.2,
                 color="#1f4e79", label="Measured (mean +/- std)")
    ax1.plot(snr_db, evm_theory_mf, "--", lw=1.2, color="#c00000",
              label="Theory, with matched-filter gain")
    ax1.plot(snr_db, evm_theory_naive, ":", lw=1.0, color="#999999",
              label="Theory, no processing gain (reference)")
    ax1.set_xlabel("SNR (dB)")
    ax1.set_ylabel("EVM (%)")
    ax1.set_yscale("log")
    ax1.set_title("EVM against SNR")
    ax1.grid(True, which="both", alpha=0.3)
    ax1.legend(fontsize=7)

    ax2.semilogy(snr_db, np.clip(ber_meas, 1e-4, None), "o-", ms=4, lw=1.2,
                 color="#1f4e79", label="Measured")
    ax2.semilogy(snr_db, np.clip(ber_theory_mf, 1e-12, None), "--", lw=1.2,
                 color="#c00000", label="Theory, with matched-filter gain")
    ax2.semilogy(snr_db, np.clip(ber_theory_naive, 1e-12, None), ":", lw=1.0,
                 color="#999999", label="Theory, no processing gain (reference)")
    ax2.set_xlabel("SNR (dB)")
    ax2.set_ylabel("BER")
    ax2.set_title("BER against SNR")
    ax2.grid(True, which="both", alpha=0.3)
    ax2.legend(fontsize=7)

    fig.tight_layout()
    fig.savefig(f"{prefix}_snr_sweep.png", dpi=300)
    plt.close(fig)


# ===========================================================================
#  MAIN RECEIVER CHAIN
# ===========================================================================

def process_capture(t, channels, args, tx_pattern):
    fs = 1.0 / np.mean(np.diff(t))
    x = channels[args.ch - 1]
    x = x - np.mean(x)

    fc = args.fc if args.fc_fixed else estimate_carrier(x, fs, args.fc)
    clock_ratio = fc / args.fc
    rs = args.rs * clock_ratio

    cutoff = args.lpf if args.lpf else 3.0 * rs
    y = passband_to_baseband(x, fs, fc, rs, cutoff)
    y = matched_filter(y, fs / rs, args.shape)

    tau0 = oerder_meyr_offset(y, fs / rs)
    tau, _ = refine_timing(y, fs, rs, tau0)

    z_raw, t_sym = extract_symbols(y, fs, rs, tau)
    if len(z_raw) < 3:
        raise ValueError(
            f"Only {len(z_raw)} symbols recovered. Use a slower timebase "
            "so the record spans more symbols."
        )

    z, df, phi = align_phase(z_raw, t_sym)
    z = normalise_power(z)

    zr, ref, bit_errs, rot, n_tied = resolve_ambiguity(z, tx_pattern)
    zr = normalise_power(zr)

    return {
        "fs": fs, "fc": fc, "rs": rs, "clock_ratio": clock_ratio,
        "tau": tau, "tau_om": tau0, "df": df, "phi": phi,
        "y": y, "z": z, "t_sym": t_sym, "x": x, "t": t,
        "z_corrected": zr, "ref": ref, "bit_errs": bit_errs,
        "n_bits": len(zr), "rotation_deg": rot * 180, "n_tied": n_tied,
    }


def run_snr_sweep(x_clean, fs, fc, rs, args, tx_pattern):
    """
    EVM-vs-SNR sweep using synthetic noise added to one clean, measured
    capture. This explores parameter space the bench data alone cannot
    reach -- it is a software experiment on real transmitter data, not a
    new physical measurement. fc and rs are held fixed at the values
    already established from the clean capture, since carrier estimation
    on a heavily noised signal is itself unreliable at low SNR -- exactly
    the failure mode the prominence warning above exists to catch.

    SNR is defined here as signal power to added-noise power over the
    full passband record. It is used as an approximation to Eb/N0 for the
    theoretical BER curve; the two are not identical in general, since
    the exact relationship depends on pulse shape, samples per symbol,
    and receiver filtering. This approximation is stated plainly rather
    than presented as an exact equivalence.
    """
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
            tau, _ = refine_timing(y, fs, rs, tau0, n_grid=41)
            z_raw, t_sym = extract_symbols(y, fs, rs, tau)
            if len(z_raw) < 3:
                continue

            z, _, _ = align_phase(z_raw, t_sym)
            z = normalise_power(z)
            zr, ref, bit_errs, _, _ = resolve_ambiguity(z, tx_pattern)
            zr = normalise_power(zr)

            evm, _ = evm_rms(zr, ref)
            evms.append(evm)
            errs_total += bit_errs
            bits_total += len(zr)

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
    sps = fs / rs

    # Two theoretical references are computed, since a naive comparison
    # against one would be misleading here:
    #
    #   evm_theory_naive treats SNR as if EVM were read from a single raw
    #   sample with no receiver processing -- it ignores the matched
    #   filter entirely.
    #
    #   evm_theory_mf accounts for the matched filter's integration gain:
    #   the integrate-and-dump filter averages over "sps" samples per
    #   symbol before a decision is made, which reduces the effective
    #   noise variance at the decision point by roughly that same factor.
    #   This is the theoretically meaningful comparison for a real
    #   receiver, and is expected to track the measured curve far more
    #   closely than the naive line -- the naive line is included only
    #   as a reference showing what a receiver with no processing gain
    #   would achieve.
    #
    # The same processing-gain correction applies to the BER curve for
    # the same physical reason: Eb/N0 at the decision point is the
    # per-sample SNR multiplied by the number of samples integrated per
    # bit.
    evm_theory_naive = 100.0 / np.sqrt(snr_lin)
    evm_theory_mf = 100.0 / np.sqrt(snr_lin * sps)
    ber_theory_naive = 0.5 * scipy_erfc(np.sqrt(snr_lin))
    ber_theory_mf = 0.5 * scipy_erfc(np.sqrt(snr_lin * sps))

    plot_snr_sweep(snr_db_arr, evm_mean, evm_std, evm_theory_naive, evm_theory_mf,
                   ber_mean, ber_theory_naive, ber_theory_mf, args.prefix)

    with open(f"{args.prefix}_snr_sweep.txt", "w", encoding="utf-8") as f:
        f.write("EVM-vs-SNR sweep, synthetic noise on a measured capture\n")
        f.write("=" * 58 + "\n")
        f.write(f"Realizations per point: {n_real}\n")
        f.write(f"Samples per symbol (matched filter integration length): {sps:.2f}\n\n")
        f.write(f"{'SNR(dB)':>10} {'EVM_meas(%)':>13} {'EVM_std(%)':>12} "
                f"{'EVM_th_naive(%)':>16} {'EVM_th_mf(%)':>13} "
                f"{'BER_meas':>12} {'BER_th_naive':>14} {'BER_th_mf':>12}\n")
        for i in range(len(snr_db_arr)):
            f.write(f"{snr_db_arr[i]:>10.1f} {evm_mean[i]:>13.3f} "
                    f"{evm_std[i]:>12.3f} {evm_theory_naive[i]:>16.3f} "
                    f"{evm_theory_mf[i]:>13.3f} {ber_mean[i]:>12.3e} "
                    f"{ber_theory_naive[i]:>14.3e} {ber_theory_mf[i]:>12.3e}\n")

    print()
    print(f"  {args.prefix}_snr_sweep.png")
    print(f"  {args.prefix}_snr_sweep.txt")


def scipy_erfc(x):
    from scipy.special import erfc
    return erfc(x)


# ===========================================================================
#  ENTRY POINT
# ===========================================================================

def main():
    p = argparse.ArgumentParser(
        description="BPSK software receiver and EVM analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument("--csv", nargs="+", required=True,
                   help="One or more scope CSV files. Several files are "
                        "pooled to improve the EVM statistics. The first "
                        "file is also used as the noise-free basis for "
                        "--snr-sweep.")
    p.add_argument("--fc", type=float, default=FC_NOMINAL_HZ,
                   help=f"Nominal carrier frequency in Hz. PLACEHOLDER "
                        f"default {FC_NOMINAL_HZ:.0f} Hz -- confirm against "
                        f"your firmware.")
    p.add_argument("--fc-fixed", action="store_true",
                   help="Trust --fc exactly instead of estimating it.")
    p.add_argument("--rs", type=float, default=RS_NOMINAL_HZ,
                   help=f"Nominal symbol rate in Bd. PLACEHOLDER default "
                        f"{RS_NOMINAL_HZ:.0f} Bd -- confirm against your "
                        f"firmware.")
    p.add_argument("--pattern", type=str, default=None,
                   help="Comma-separated transmitted bit pattern, e.g. "
                        "0,1,0,1,0,1,0,1. REQUIRED for a trustworthy "
                        "result -- if omitted, a placeholder pattern is "
                        "used and the result is not to be trusted.")
    p.add_argument("--shape", choices=["nrz", "rc"], default="nrz",
                   help="Transmit pulse shape, must match the firmware.")
    p.add_argument("--lpf", type=float, default=None,
                   help="Low-pass cutoff in Hz (default 3 x Rs).")
    p.add_argument("--ch", type=int, default=1,
                   help="Which channel holds the signal.")
    p.add_argument("--prefix", default="bpsk",
                   help="Output filename prefix.")
    p.add_argument("--exclude-ambiguous", action="store_true",
                   help="Drop any capture whose rotation/shift alignment "
                        "is not uniquely determined (tie count > 1) before "
                        "pooling, rather than including it anyway. "
                        "Recommended whenever the per-capture symbol count "
                        "is small relative to the pattern length.")
    p.add_argument("--no-plots", action="store_true")

    p.add_argument("--snr-sweep", action="store_true",
                   help="Run an EVM-vs-SNR sweep using the first --csv "
                        "file as the noise-free basis, instead of (or in "
                        "addition to) the direct pooled analysis.")
    p.add_argument("--snr-min", type=float, default=-5.0)
    p.add_argument("--snr-max", type=float, default=20.0)
    p.add_argument("--snr-step", type=float, default=1.0)
    p.add_argument("--snr-realizations", type=int, default=30,
                   help="Independent noise draws averaged per SNR point.")
    p.add_argument("--snr-seed", type=int, default=0)

    args = p.parse_args()

    if args.pattern is None:
        print("=" * 68)
        print("WARNING: --pattern was not supplied. Using a PLACEHOLDER "
              "bit pattern.")
        print("This result is NOT trustworthy until you supply the actual "
              "transmitted pattern from your firmware, e.g.:")
        print("  --pattern 0,1,0,1,0,1,0,1")
        print("=" * 68)
        tx_pattern = TX_PATTERN_DEFAULT
    else:
        tx_pattern = np.array([int(b) for b in args.pattern.split(",")], dtype=int)
        if not np.all((tx_pattern == 0) | (tx_pattern == 1)):
            print("ERROR: --pattern must contain only 0s and 1s.")
            sys.exit(1)

    print("=" * 68)
    print("BPSK SOFTWARE RECEIVER AND EVM ANALYSIS")
    print("=" * 68)
    print(f"Pulse shape   : {args.shape}")
    print(f"Nominal fc    : {args.fc:.1f} Hz")
    print(f"Nominal Rs    : {args.rs:.1f} Bd")
    print(f"TX pattern    : {list(tx_pattern)}")
    print()
    print("Loading captures")

    results = []
    for path in args.csv:
        t, chans = load_tek_csv(path)
        try:
            r = process_capture(t, chans, args, tx_pattern)
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
        tie_note = f", ties={r['n_tied']}" if r['n_tied'] > 1 else ""
        print(f"  {r['file']}: rotation={r['rotation_deg']:>3.0f} deg, "
              f"symbols={len(r['z_corrected'])}, "
              f"bit errors={r['bit_errs']} of {r['n_bits']}{tie_note}")

    if args.exclude_ambiguous:
        kept = [r for r in results if r["n_tied"] <= 1]
        dropped = [r for r in results if r["n_tied"] > 1]
        if dropped:
            print()
            print(f"--exclude-ambiguous: dropping {len(dropped)} of "
                  f"{len(results)} captures with an undetermined alignment:")
            for r in dropped:
                print(f"    {r['file']} (ties={r['n_tied']})")
        if not kept:
            print()
            print("ERROR: --exclude-ambiguous left no usable captures. "
                  "Every capture's alignment was ambiguous -- the record "
                  "is too short relative to the pattern length for this "
                  "dataset to support a reliable result. A longer capture "
                  "(slower timebase) is needed.")
            sys.exit(1)
        results = kept

    z_all = np.concatenate([r["z_corrected"] for r in results])
    ref_all = np.concatenate([r["ref"] for r in results])
    z_all = normalise_power(z_all)

    evm, err_vec = evm_rms(z_all, ref_all)
    evm_peak = 100.0 * np.max(np.abs(err_vec)) / np.sqrt(np.mean(np.abs(IDEAL[ref_all]) ** 2))
    mer_db = -20.0 * np.log10(evm / 100.0)

    ph_rms, ph_pk = phase_error_stats(z_all, ref_all)
    resid_q = residual_quadrature(z_all)

    bit_errs = int(sum(r["bit_errs"] for r in results))
    n_bits = int(sum(r["n_bits"] for r in results))

    r0 = results[0]
    zr = z_all

    print()
    print("-" * 68)
    print("RECOVERY")
    print("-" * 68)
    ppm = (r0["clock_ratio"] - 1.0) * 1e6
    print(f"Carrier estimate      : {r0['fc']:.4f} Hz ({ppm:+.0f} ppm vs nominal)")
    print(f"Symbol rate (scaled)  : {r0['rs']:.4f} Bd")
    print(f"Scope sample rate     : {r0['fs'] / 1e3:.3f} kSa/s")
    print(f"Samples per symbol    : {r0['fs'] / r0['rs']:.2f}")
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
    print(f"Residual quadrature   : {resid_q:.3f} % (diagnostic; should be small)")
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
        plot_spectrum(r0["x"], r0["fs"], args.prefix, fc=r0["fc"], rs=r0["rs"])
        plot_iq_tracks(r0["y"], r0["fs"], r0["t_sym"], r0["z"], args.prefix)
        plot_constellation(zr, ref_all, evm, args.prefix)
        plot_eye(r0["y"], r0["fs"], r0["rs"], r0["tau"], args.prefix)
        for suffix in ("waveform", "spectrum", "iq_tracks", "constellation", "eye"):
            f = f"{args.prefix}_{suffix}.png"
            if os.path.isfile(f):
                print(f"  {f}")

    with open(f"{args.prefix}_results.txt", "w", encoding="utf-8") as f:
        f.write("BPSK software receiver and EVM analysis\n")
        f.write("=" * 44 + "\n")
        f.write(f"Pulse shape           : {args.shape}\n")
        f.write(f"Files                 : {', '.join(args.csv)}\n")
        f.write(f"TX pattern            : {list(tx_pattern)}\n")
        f.write(f"Carrier estimate      : {r0['fc']:.4f} Hz\n")
        f.write(f"Symbol rate           : {r0['rs']:.4f} Bd\n")
        f.write(f"Scope sample rate     : {r0['fs']:.1f} Sa/s\n")
        f.write(f"Symbols analysed      : {len(zr)}\n")
        f.write(f"EVM (RMS)             : {evm:.3f} %\n")
        f.write(f"EVM (peak)            : {evm_peak:.3f} %\n")
        f.write(f"MER                   : {mer_db:.2f} dB\n")
        f.write(f"Phase error, RMS      : {ph_rms:.3f} deg\n")
        f.write(f"Bit errors            : {bit_errs} of {n_bits}\n")
    print(f"  {args.prefix}_results.txt")

    if args.snr_sweep:
        run_snr_sweep(r0["x"], r0["fs"], r0["fc"], r0["rs"], args, tx_pattern)


if __name__ == "__main__":
    main()
