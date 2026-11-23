"""
make_figures.py
================
Generates every Chapter 7 figure for the phase deviation and phase
accuracy metric that does not already exist, from CSV captures already
taken. No new bench work is required.

Figures 7.42 (transition detail) and 7.43 (zero-crossing verification)
are NOT regenerated here -- they are produced by m3_transition_detail.py
and m3_simple_demo.py respectively and already exist.

Figures 7.30 and 7.40 are plotted from the exported CSV records rather
than photographed from the instrument screen, which gives cleaner
reproduction in print and avoids camera glare.

Usage:
    1. Edit the CONFIG block below so the filenames match your files.
    2. python make_figures.py

Each figure is generated independently. If a file is missing the script
reports it and continues, so a partial run still produces everything it
can.

Requires in the same folder: measure_fc_ref.py, check_m1_point.py,
check_m3_point.py
"""

import os
import traceback
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import hilbert, butter, filtfilt

from measure_fc_ref import load_tbs_csv
from check_m1_point import (
    extract_phase, measure_beta, fit_delay_offset, wrap,
    TRIM_CYCLES, BW_HZ, FM_HZ,
)
from check_m3_point import (
    find_sync_rising_edge, bpsk_ideal_phase, steady_state_mask,
    RS, GUARD, bits,
)

# ==========================================================================
#                                CONFIG
# ==========================================================================
FC_HAT = 5000.422114          # from fc_ref.json (R1 zero-crossing estimate)

OUTDIR = "figures"

# --- reference captures ---
REF_CARRIER = "ref_carrier_kp0.csv"      # R1
REF_NOISE = "R3.CSV"                     # R3

# --- the kp = 1.0 rad capture, used for figures 7.30 - 7.34 ---
KP1_CH1 = "m1_kp1p00_r1.csv"
KP1_CH2 = None       # set to a filename if a CH2 record exists for this
                     # capture; leave None if only CH1 was saved

# --- M1 sweep ---
M1_MANIFEST = "m1_manifest.txt"

# --- M2 fine sweep ---
M2_FINE_KPS = [2.30, 2.32, 2.34, 2.36, 2.38, 2.40,
               2.42, 2.44, 2.46, 2.48, 2.50]
M2_FINE_PATTERN = "m2_j0_fine_{kp:.2f}.csv"
M2_NULL_FILE = "m2_j0_fine_2.42.csv"     # capture nearest the fitted null,
                                          # used for the spectrum comparison

# --- M3 ---
M3_CH1 = "m3_ch1_r1.csv"
M3_CH2 = "m3_ch2_r1.csv"
M3_BW = 2500.0

J0_NULL_THEORY = 2.40483
FULL_SCALE_RAD = 4.0
# ==========================================================================


def save(fig, name):
    os.makedirs(OUTDIR, exist_ok=True)
    path = os.path.join(OUTDIR, name)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"  saved {path}")


def carrier_amplitude(t, x, f_carrier):
    xa = x - np.mean(x)
    N = len(xa)
    dt = (t[-1] - t[0]) / (N - 1)
    w = np.blackman(N)
    X = np.abs(np.fft.rfft(xa * w))
    f = np.fft.rfftfreq(N, dt)
    k = int(np.argmin(np.abs(f - f_carrier)))
    return float(X[max(k - 2, 0):k + 3].max())


def spectrum_db(t, x, ref_level):
    xa = x - np.mean(x)
    N = len(xa)
    dt = (t[-1] - t[0]) / (N - 1)
    w = np.blackman(N)
    X = np.abs(np.fft.rfft(xa * w))
    f = np.fft.rfftfreq(N, dt)
    with np.errstate(divide="ignore"):
        db = 20 * np.log10(np.maximum(X, 1e-12) / ref_level)
    return f, db


def read_manifest(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            kp_str, fname = line.split(None, 1)
            rows.append((float(kp_str), fname.strip()))
    return rows


# ==========================================================================
# 7.30 -- PM output waveform (and sync marker if available)
# ==========================================================================
def fig_7_30():
    t, x = load_tbs_csv(KP1_CH1, channel=1)

    if KP1_CH2:
        t2, x2 = load_tbs_csv(KP1_CH2, channel=1)
        fig, ax = plt.subplots(2, 1, figsize=(9, 5), sharex=True,
                               gridspec_kw={"height_ratios": [3, 1]})
        ax[0].plot(t * 1e3, x, color="tab:blue", lw=0.9)
        ax[0].set_ylabel("LPF output (V)")
        ax[0].grid(True, alpha=0.3)
        ax[1].plot(t2 * 1e3, x2, color="tab:orange", lw=0.9)
        ax[1].set_ylabel("Sync (V)")
        ax[1].set_xlabel("Time (ms)")
        ax[1].grid(True, alpha=0.3)
        ax[0].set_title("Phase-modulated output at $k_p$ = 1.0 rad, "
                        "with synchronisation marker")
    else:
        fig, ax0 = plt.subplots(figsize=(9, 3.5))
        ax0.plot(t * 1e3, x, color="tab:blue", lw=0.9)
        ax0.set_ylabel("LPF output (V)")
        ax0.set_xlabel("Time (ms)")
        ax0.grid(True, alpha=0.3)
        ax0.set_title("Phase-modulated output at $k_p$ = 1.0 rad")
        print("  note: no CH2 file configured -- sync marker panel omitted. "
              "Adjust the figure caption accordingly.")

    save(fig, "fig_7_30_pm_waveform.png")


# ==========================================================================
# 7.31 -- envelope
# ==========================================================================
def fig_7_31():
    t, x = load_tbs_csv(KP1_CH1, channel=1)
    fs = 1.0 / ((t[-1] - t[0]) / (len(t) - 1))
    xa = x - np.mean(x)
    zb = hilbert(xa) * np.exp(-1j * 2 * np.pi * FC_HAT * t)
    b, a = butter(4, BW_HZ / (fs / 2.0))
    zb = filtfilt(b, a, zb.real) + 1j * filtfilt(b, a, zb.imag)
    n = int(TRIM_CYCLES * fs / FC_HAT)
    tt, env = t[n:-n], np.abs(zb)[n:-n]
    ripple = (env.max() - env.min()) / env.mean() * 100

    fig, ax = plt.subplots(figsize=(9, 3.5))
    ax.plot(tt * 1e3, env, color="tab:green", lw=1.2)
    ax.axhline(env.mean(), color="gray", ls=":", lw=1, label="Mean")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Envelope |$z_b(t)$| (V)")
    ax.set_title(f"Instantaneous envelope at $k_p$ = 1.0 rad "
                 f"(ripple {ripple:.2f} %)")
    ax.set_ylim(0, env.max() * 1.25)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    save(fig, "fig_7_31_envelope.png")


# ==========================================================================
# 7.32 / 7.33 / 7.34 -- trajectory, error, histogram
# ==========================================================================
def fig_7_32_33_34():
    kp = 1.0
    t, x = load_tbs_csv(KP1_CH1, channel=1)
    t_p, theta, _ = extract_phase(t, x, FC_HAT, BW_HZ, TRIM_CYCLES)
    ideal = lambda tt: kp * np.sin(2 * np.pi * FM_HZ * tt)
    tau, offset, err = fit_delay_offset(t_p, theta, ideal)
    err_deg = np.degrees(err)
    rms = float(np.sqrt(np.mean(err_deg ** 2)))
    peak = float(np.max(np.abs(err_deg)))

    fig, ax = plt.subplots(figsize=(9, 3.5))
    ax.plot(t_p * 1e3, theta, color="tab:blue", lw=1.2, label="Measured")
    ax.plot(t_p * 1e3, ideal(t_p - tau) + offset, "--", color="tab:red",
            lw=1.0, label="Ideal (delay and offset removed)")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Phase deviation (rad)")
    ax.set_title("Measured and ideal phase trajectories, "
                 "$k_p$ = 1.0 rad")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    save(fig, "fig_7_32_trajectory.png")

    fig, ax = plt.subplots(figsize=(9, 3))
    ax.plot(t_p * 1e3, err_deg, color="tab:purple", lw=0.9)
    ax.axhline(0, color="gray", ls=":", lw=1)
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Phase error (deg)")
    ax.set_title(f"Phase error against time, $k_p$ = 1.0 rad "
                 f"(RMS {rms:.3f}\u00b0, peak {peak:.3f}\u00b0, "
                 f"post-detection bandwidth {BW_HZ/1e3:.1f} kHz)")
    ax.grid(True, alpha=0.3)
    save(fig, "fig_7_33_phase_error.png")

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(err_deg, bins=60, color="tab:blue", edgecolor="black", lw=0.4)
    ax.set_xlabel("Phase error (deg)")
    ax.set_ylabel("Number of samples")
    ax.set_title(f"Distribution of phase error, $k_p$ = 1.0 rad "
                 f"(RMS {rms:.3f}\u00b0)")
    ax.grid(True, alpha=0.3)
    save(fig, "fig_7_34_histogram.png")


# ==========================================================================
# 7.35 / 7.36 / 7.37 -- M1 sweep, residuals, gain error
# ==========================================================================
def fig_7_35_36_37():
    rows = read_manifest(M1_MANIFEST)
    kps, betas = [], []
    for kp, fname in rows:
        t, x = load_tbs_csv(fname, channel=1)
        t_p, theta, _ = extract_phase(t, x, FC_HAT, BW_HZ, TRIM_CYCLES)
        kps.append(kp)
        betas.append(measure_beta(t_p, theta, FM_HZ))
    kps = np.array(kps)
    betas = np.array(betas)

    A = np.vstack([kps, np.ones_like(kps)]).T
    coef, *_ = np.linalg.lstsq(A, betas, rcond=None)
    slope, intercept = coef
    pred = A @ coef
    resid = betas - pred
    r2 = 1 - np.sum(resid ** 2) / np.sum((betas - betas.mean()) ** 2)
    lin_pct = np.max(np.abs(resid)) / FULL_SCALE_RAD * 100

    # 7.35 -- the key figure
    fig, ax = plt.subplots(figsize=(7, 5.5))
    xs = np.linspace(0, max(kps) * 1.05, 100)
    ax.plot(xs, slope * xs + intercept, "-", color="tab:red", lw=1.2,
            label=f"Fit: slope {slope:.4f}, $R^2$ = {r2:.6f}")
    ax.plot(xs, xs, ":", color="gray", lw=1.0, label="Ideal (slope 1)")
    ax.plot(kps, betas, "o", color="tab:blue", ms=5, label="Measured")
    ax.set_xlabel("Commanded phase deviation $k_p$ (rad)")
    ax.set_ylabel(r"Measured phase deviation $\hat{\beta}$ (rad)")
    ax.set_title(f"Commanded against measured phase deviation "
                 f"({len(rows)} captures)")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(True, alpha=0.3)
    ax.set_aspect("equal", adjustable="box")
    save(fig, "fig_7_35_m1_sweep.png")

    # 7.36 -- residuals
    fig, ax = plt.subplots(figsize=(9, 3.5))
    ax.plot(kps, resid * 1e3, "o", color="tab:blue", ms=5)
    ax.axhline(0, color="gray", ls=":", lw=1)
    ax.set_xlabel("Commanded phase deviation $k_p$ (rad)")
    ax.set_ylabel("Residual (mrad)")
    ax.set_title(f"Residuals of the linear fit "
                 f"(maximum deviation {lin_pct:.2f} % of full scale)")
    ax.grid(True, alpha=0.3)
    save(fig, "fig_7_36_m1_residuals.png")

    # 7.37 -- gain error against kp, with spread across repeats
    uniq = sorted(set(kps))
    means, sds = [], []
    for u in uniq:
        g = (betas[kps == u] / u - 1) * 100
        means.append(g.mean())
        sds.append(g.std())
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.errorbar(uniq, means, yerr=sds, fmt="o-", color="tab:blue",
                ms=5, lw=1.0, capsize=3)
    ax.axhline(0, color="gray", ls=":", lw=1)
    ax.set_xlabel("Commanded phase deviation $k_p$ (rad)")
    ax.set_ylabel("Gain error (%)")
    ax.set_title("Gain error against commanded phase deviation "
                 "(error bars: spread across repeats)")
    ax.grid(True, alpha=0.3)
    save(fig, "fig_7_37_gain_error.png")


# ==========================================================================
# 7.38 -- M2 fine sweep with parabolic fit
# ==========================================================================
def fig_7_38():
    t_ref, x_ref = load_tbs_csv(REF_CARRIER, channel=1)
    amp_ref = carrier_amplitude(t_ref, x_ref, FC_HAT)

    kps, dbs = [], []
    for kp in M2_FINE_KPS:
        t, x = load_tbs_csv(M2_FINE_PATTERN.format(kp=kp), channel=1)
        a = carrier_amplitude(t, x, FC_HAT)
        kps.append(kp)
        dbs.append(20 * np.log10(a / amp_ref))
    kps = np.array(kps)
    dbs = np.array(dbs)

    # parabola fitted to log magnitude, as in m2_j0_fine_fit.py
    log_amp = dbs / 20.0
    pa, pb, pc = np.polyfit(kps, log_amp, 2)
    kp_null = -pb / (2 * pa)
    gain_err = (kp_null - J0_NULL_THEORY) / J0_NULL_THEORY * 100
    xs = np.linspace(kps.min(), kps.max(), 300)
    fit_db = np.polyval([pa, pb, pc], xs) * 20

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(xs, fit_db, "-", color="tab:red", lw=1.2, label="Parabolic fit")
    ax.plot(kps, dbs, "o", color="tab:blue", ms=5, label="Measured")
    ax.axvline(kp_null, color="tab:green", ls="--", lw=1.0,
               label=f"Fitted null: {kp_null:.4f} rad")
    ax.axvline(J0_NULL_THEORY, color="gray", ls=":", lw=1.0,
               label=f"Theory: {J0_NULL_THEORY:.5f} rad")
    ax.set_xlabel("Commanded phase sensitivity $k_p$ (rad)")
    ax.set_ylabel("Carrier amplitude relative to unmodulated (dB)")
    ax.set_title(f"First Bessel carrier null "
                 f"(gain error {gain_err:+.2f} %)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    save(fig, "fig_7_38_m2_fine.png")


# ==========================================================================
# 7.39 -- spectra, unmodulated against the null
# ==========================================================================
def fig_7_39():
    t_ref, x_ref = load_tbs_csv(REF_CARRIER, channel=1)
    amp_ref = carrier_amplitude(t_ref, x_ref, FC_HAT)
    f_ref, db_ref = spectrum_db(t_ref, x_ref, amp_ref)

    t_n, x_n = load_tbs_csv(M2_NULL_FILE, channel=1)
    f_n, db_n = spectrum_db(t_n, x_n, amp_ref)

    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.plot(f_ref / 1e3, db_ref, color="gray", lw=1.0,
            label="Unmodulated ($k_p$ = 0)")
    ax.plot(f_n / 1e3, db_n, color="tab:blue", lw=1.0,
            label="At the carrier null")
    ax.axvline(FC_HAT / 1e3, color="tab:red", ls=":", lw=1.0,
               label="$f_c$")
    ax.set_xlim(0, 15)
    ax.set_ylim(-70, 5)
    ax.set_xlabel("Frequency (kHz)")
    ax.set_ylabel("Magnitude relative to unmodulated carrier (dB)")
    ax.set_title("Spectrum at the first Bessel carrier null "
                 "(Blackman window)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    save(fig, "fig_7_39_spectra.png")


# ==========================================================================
# 7.40 -- BPSK waveform with sync marker
# ==========================================================================
def fig_7_40():
    t1, x1 = load_tbs_csv(M3_CH1, channel=1)
    t2, x2 = load_tbs_csv(M3_CH2, channel=1)

    fig, ax = plt.subplots(2, 1, figsize=(9, 5), sharex=True,
                           gridspec_kw={"height_ratios": [3, 1]})
    ax[0].plot(t1 * 1e3, x1, color="tab:blue", lw=0.7)
    ax[0].set_ylabel("LPF output (V)")
    ax[0].set_title("BPSK output showing commanded phase reversals, "
                    "with synchronisation marker")
    ax[0].grid(True, alpha=0.3)
    ax[1].plot(t2 * 1e3, x2, color="tab:orange", lw=0.9)
    ax[1].set_ylabel("Sync (V)")
    ax[1].set_xlabel("Time (ms)")
    ax[1].grid(True, alpha=0.3)
    save(fig, "fig_7_40_bpsk_waveform.png")


# ==========================================================================
# 7.41 -- BPSK phase trajectory with guard intervals
# ==========================================================================
def fig_7_41():
    t1, x1 = load_tbs_csv(M3_CH1, channel=1)
    t2, x2 = load_tbs_csv(M3_CH2, channel=1)
    t_sync = find_sync_rising_edge(t2, x2)

    t_p, theta, _ = extract_phase(t1, x1, FC_HAT, M3_BW, TRIM_CYCLES)
    ideal = bpsk_ideal_phase(t_p, t_sync, RS, bits)
    mask = steady_state_mask(t_p, t_sync, RS, GUARD)

    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(t_p * 1e3, np.degrees(theta), color="tab:blue", lw=1.0,
            label="Measured")
    ax.plot(t_p * 1e3, np.degrees(ideal), "--", color="tab:red", lw=1.0,
            label="Ideal")

    # shade the guard intervals
    excluded = ~mask
    edges = np.diff(excluded.astype(int))
    starts = np.where(edges == 1)[0]
    ends = np.where(edges == -1)[0]
    if excluded[0]:
        starts = np.insert(starts, 0, 0)
    if excluded[-1]:
        ends = np.append(ends, len(excluded) - 1)
    for i, (s, e) in enumerate(zip(starts, ends)):
        ax.axvspan(t_p[s] * 1e3, t_p[e] * 1e3, color="gray", alpha=0.2,
                   label="Guard interval" if i == 0 else None)

    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Phase (deg)")
    ax.set_title(f"BPSK phase trajectory with "
                 f"\u00b1{GUARD*1e6:.0f} \u00b5s guard intervals "
                 f"(post-detection bandwidth {M3_BW/1e3:.1f} kHz)")
    ax.legend(fontsize=8, loc="center right")
    ax.grid(True, alpha=0.3)
    save(fig, "fig_7_41_bpsk_trajectory.png")


# ==========================================================================
# 7.44 -- noise record and spectrum
# ==========================================================================
def fig_7_44():
    t, x = load_tbs_csv(REF_NOISE, channel=1)
    dt = (t[-1] - t[0]) / (len(t) - 1)
    xa = x - np.mean(x)
    rms_mv = float(np.std(xa)) * 1e3

    N = len(xa)
    w = np.blackman(N)
    X = np.abs(np.fft.rfft(xa * w))
    f = np.fft.rfftfreq(N, dt)
    with np.errstate(divide="ignore"):
        db = 20 * np.log10(np.maximum(X, 1e-12) / X.max())

    fig, ax = plt.subplots(2, 1, figsize=(9, 6))
    ax[0].plot(t * 1e3, xa * 1e3, color="tab:blue", lw=0.6)
    ax[0].set_xlabel("Time (ms)")
    ax[0].set_ylabel("Noise (mV)")
    ax[0].set_title(f"Noise floor with the DAC held at mid-scale "
                    f"(RMS {rms_mv:.2f} mV)")
    ax[0].grid(True, alpha=0.3)

    ax[1].plot(f / 1e3, db, color="tab:purple", lw=0.8)
    ax[1].set_xlim(0, 20)
    ax[1].set_xlabel("Frequency (kHz)")
    ax[1].set_ylabel("Relative magnitude (dB)")
    ax[1].set_title("Noise spectrum")
    ax[1].grid(True, alpha=0.3)
    save(fig, "fig_7_44_noise.png")


# ==========================================================================
FIGURES = [
    ("7.30  PM waveform", fig_7_30),
    ("7.31  envelope", fig_7_31),
    ("7.32/33/34  trajectory, error, histogram", fig_7_32_33_34),
    ("7.35/36/37  M1 sweep, residuals, gain error", fig_7_35_36_37),
    ("7.38  M2 fine sweep", fig_7_38),
    ("7.39  spectra", fig_7_39),
    ("7.40  BPSK waveform", fig_7_40),
    ("7.41  BPSK trajectory", fig_7_41),
    ("7.44  noise floor", fig_7_44),
]


def main():
    print(f"Output directory: {OUTDIR}/")
    print(f"Carrier frequency: {FC_HAT:.6f} Hz\n")
    ok, failed = 0, []
    for label, fn in FIGURES:
        print(f"[{label}]")
        try:
            fn()
            ok += 1
        except FileNotFoundError as e:
            print(f"  SKIPPED -- file not found: {e.filename}")
            failed.append(label)
        except Exception:
            print("  FAILED:")
            traceback.print_exc()
            failed.append(label)
        print()

    print(f"=== {ok} of {len(FIGURES)} figure groups generated ===")
    if failed:
        print("Not generated:")
        for f in failed:
            print(f"  - {f}")
        print("\nCheck the filenames in the CONFIG block at the top of "
              "this script against your actual files.")
    print("\nFigures 7.42 and 7.43 are produced by m3_transition_detail.py "
          "and m3_simple_demo.py and are not regenerated here.")


if __name__ == "__main__":
    main()
