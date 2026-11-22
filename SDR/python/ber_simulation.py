#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ber_simulation.py
==================

BPSK bit-error-rate simulation against theory, over an AWGN channel.

Project : Implementation and Performance Analysis of a Baseband Phase
          Modulation System with STM32 for Low-Power Wireless Communication

Purpose
-------
This is a deliberately idealized, perfectly-synchronized Monte Carlo
simulation -- the theoretical ceiling that the real software receiver
(bpsk_analysis.py) is compared against. It exists to answer a question
the real hardware captures cannot: what is the true BER at a given
Eb/N0, with enough bits simulated to make the number statistically
trustworthy? A handful of real captures can show a *trend* (see the
EVM-vs-SNR sweep in bpsk_analysis.py); this script is what supplies an
actual, resolvable BER curve.

No timing recovery, no carrier recovery, no matched filtering delay --
symbols are decided directly, since the point of this simulation is the
theoretical limit, not a receiver implementation. bpsk_analysis.py is
where receiver imperfections are characterised; this script deliberately
excludes them.

Method
------
For each Eb/N0 point:
  1. Generate random bits, map to +1 / -1 (unit energy per bit, matching
     the project's BPSK convention: bit 0 -> +1, bit 1 -> -1).
  2. Add real AWGN with variance N0/2, derived from the target Eb/N0
     (Eb = 1 by construction, so N0 = 1 / Eb_N0_linear).
  3. Decide each symbol by the sign of the received sample (optimal for
     BPSK in AWGN).
  4. Count bit errors.
  5. Repeat in batches until EITHER a minimum number of bit errors has
     accumulated (for statistical reliability) OR a maximum bit ceiling
     is reached (so low-BER points near high Eb/N0 do not run forever).

Validation
----------
At Eb/N0 = 0 dB, theoretical BPSK BER is approximately 0.0786. The
script checks the simulated value against this automatically and prints
a pass/fail line -- this is the sanity check that must pass before any
other result from this script is trusted.

Usage
-----
  python ber_simulation.py
  python ber_simulation.py --ebno-min 0 --ebno-max 12 --ebno-step 0.5 \\
      --min-errors 100 --max-bits 1e8 --prefix ber_bpsk
"""

import argparse

import numpy as np
from scipy.special import erfc

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ===========================================================================
#  THEORY
# ===========================================================================

def q_function(x):
    """Gaussian Q-function via the complementary error function."""
    return 0.5 * erfc(x / np.sqrt(2.0))


def bpsk_theory_ber(ebno_db):
    """Theoretical BPSK BER in AWGN: BER = Q(sqrt(2 * Eb/N0))."""
    ebno_lin = 10.0 ** (np.asarray(ebno_db) / 10.0)
    return q_function(np.sqrt(2.0 * ebno_lin))


# ===========================================================================
#  SIMULATION
# ===========================================================================

def simulate_bpsk_ber(ebno_db, min_errors=100, max_bits=int(1e8),
                       batch_size=200_000, rng=None):
    """
    Monte Carlo BPSK BER at one Eb/N0 point, with an error-count stopping
    rule: accumulate bits in batches until either min_errors bit errors
    have been observed, or max_bits total bits have been simulated,
    whichever comes first.

    Eb is fixed at 1 by construction (unit-energy symbols), so the noise
    variance is derived directly from the target Eb/N0:
        N0 = Eb / (Eb/N0)_linear = 1 / (Eb/N0)_linear
        noise variance (real, one-sided) = N0 / 2

    Returns a dict with the bit count, error count, empirical BER, and a
    Wilson score 95% confidence interval on the BER estimate.
    """
    if rng is None:
        rng = np.random.default_rng()

    ebno_lin = 10.0 ** (ebno_db / 10.0)
    noise_std = np.sqrt(1.0 / (2.0 * ebno_lin))

    bits_total = 0
    errors_total = 0

    while errors_total < min_errors and bits_total < max_bits:
        n = min(batch_size, max_bits - bits_total)

        bits = rng.integers(0, 2, size=n)
        tx = np.where(bits == 0, 1.0, -1.0)          # bit 0 -> +1, bit 1 -> -1
        noise = rng.normal(0.0, noise_std, size=n)
        rx = tx + noise

        dec_bits = np.where(rx >= 0.0, 0, 1)
        errors_total += int(np.sum(dec_bits != bits))
        bits_total += n

    ber = errors_total / bits_total if bits_total > 0 else np.nan
    lo, hi = wilson_ci(errors_total, bits_total)

    return {
        "ebno_db": ebno_db, "bits": bits_total, "errors": errors_total,
        "ber": ber, "ci_lo": lo, "ci_hi": hi,
        "hit_error_target": errors_total >= min_errors,
    }


def wilson_ci(k, n, z=1.96):
    """
    Wilson score 95% confidence interval for a binomial proportion
    (here, bit-error rate from k errors in n trials). Used in preference
    to a plain normal approximation because it stays well-behaved even
    when k is small or the true rate is very low -- exactly the regime
    most points on a BER curve fall into.
    """
    if n == 0:
        return 0.0, 0.0
    phat = k / n
    denom = 1.0 + z ** 2 / n
    center = phat + z ** 2 / (2.0 * n)
    adj = z * np.sqrt(phat * (1.0 - phat) / n + z ** 2 / (4.0 * n ** 2))
    return max(0.0, (center - adj) / denom), min(1.0, (center + adj) / denom)


# ===========================================================================
#  PLOTTING
# ===========================================================================

def plot_ber_curve(results, prefix):
    ebno = np.array([r["ebno_db"] for r in results])
    ber = np.array([r["ber"] for r in results])
    ci_lo = np.array([r["ci_lo"] for r in results])
    ci_hi = np.array([r["ci_hi"] for r in results])

    ebno_fine = np.linspace(ebno.min(), ebno.max(), 400)
    theory = bpsk_theory_ber(ebno_fine)

    yerr_lo = np.clip(ber - ci_lo, 0, None)
    yerr_hi = np.clip(ci_hi - ber, 0, None)

    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    ax.semilogy(ebno_fine, theory, "-", lw=1.4, color="#c00000",
                label="Theory, Q(sqrt(2 Eb/N0))")
    ax.errorbar(ebno, np.clip(ber, 1e-12, None), yerr=[yerr_lo, yerr_hi],
                fmt="o", ms=5, lw=1.0, capsize=3, color="#1f4e79",
                label="Simulated (95% Wilson CI)")

    ax.set_xlabel("Eb/N0 (dB)")
    ax.set_ylabel("Bit error rate")
    ax.set_title("BPSK BER against Eb/N0 (AWGN, idealised receiver)")
    ax.set_ylim(bottom=max(1e-9, np.nanmin(ber[ber > 0]) * 0.3) if np.any(ber > 0) else 1e-9)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(f"{prefix}_curve.png", dpi=300)
    plt.close(fig)


# ===========================================================================
#  MAIN
# ===========================================================================

def main():
    p = argparse.ArgumentParser(description="BPSK BER-vs-Eb/N0 Monte Carlo simulation")
    p.add_argument("--ebno-min", type=float, default=0.0)
    p.add_argument("--ebno-max", type=float, default=12.0)
    p.add_argument("--ebno-step", type=float, default=1.0)
    p.add_argument("--min-errors", type=int, default=100,
                   help="Stop a point once this many bit errors have "
                        "accumulated (statistical reliability target).")
    p.add_argument("--max-bits", type=float, default=1e8,
                   help="Bit ceiling per point, so low-BER points at high "
                        "Eb/N0 cannot run indefinitely.")
    p.add_argument("--batch-size", type=int, default=200_000)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--prefix", default="ber_bpsk")
    args = p.parse_args()

    rng = np.random.default_rng(args.seed)
    ebno_range = np.arange(args.ebno_min, args.ebno_max + 1e-9, args.ebno_step)

    print("=" * 68)
    print("BPSK BER-VS-Eb/N0 MONTE CARLO SIMULATION")
    print("=" * 68)
    print(f"Eb/N0 range      : {args.ebno_min:.1f} to {args.ebno_max:.1f} dB, "
          f"step {args.ebno_step:.1f} dB")
    print(f"Stopping rule    : {args.min_errors} errors, or "
          f"{args.max_bits:.0e} bits, whichever first")
    print()

    # ---- Validation check (Part E1): Eb/N0 = 0 dB should give BER ~ 0.0786 ----
    print("-" * 68)
    print("VALIDATION CHECK -- must pass before trusting any result below")
    print("-" * 68)
    val = simulate_bpsk_ber(0.0, min_errors=max(args.min_errors, 200),
                             max_bits=int(args.max_bits), rng=rng)
    theory_val = float(bpsk_theory_ber(0.0))
    rel_err = abs(val["ber"] - theory_val) / theory_val
    passed = rel_err < 0.15   # generous given finite-sample noise at ~200 errors
    print(f"  Eb/N0 = 0 dB: simulated BER = {val['ber']:.4f} "
          f"({val['errors']} errors / {val['bits']} bits), "
          f"theory = {theory_val:.4f}")
    print(f"  Relative difference: {rel_err*100:.1f}%  ->  "
          f"{'PASS' if passed else 'FAIL -- DO NOT TRUST RESULTS BELOW'}")
    if not passed:
        print("  Stopping: fix the simulation before proceeding.")
        return
    print()

    # ---- Full sweep ----
    print("-" * 68)
    print("SWEEP")
    print("-" * 68)
    print(f"{'Eb/N0(dB)':>10} {'BER (sim)':>12} {'errors':>8} {'bits':>12} "
          f"{'BER (theory)':>13} {'ratio':>8}")

    results = []
    for ebno_db in ebno_range:
        r = simulate_bpsk_ber(ebno_db, min_errors=args.min_errors,
                              max_bits=int(args.max_bits),
                              batch_size=args.batch_size, rng=rng)
        results.append(r)
        theory = float(bpsk_theory_ber(ebno_db))
        ratio = r["ber"] / theory if theory > 0 else np.nan
        flag = "" if r["hit_error_target"] else "  (bit ceiling hit)"
        print(f"{ebno_db:>10.1f} {r['ber']:>12.3e} {r['errors']:>8d} "
              f"{r['bits']:>12d} {theory:>13.3e} {ratio:>8.3f}{flag}")

    plot_ber_curve(results, args.prefix)

    with open(f"{args.prefix}_results.txt", "w", encoding="utf-8") as f:
        f.write("BPSK BER-vs-Eb/N0 Monte Carlo simulation\n")
        f.write("=" * 44 + "\n")
        f.write(f"Stopping rule: {args.min_errors} errors or "
                f"{args.max_bits:.0e} bits\n")
        f.write(f"Validation at 0 dB: simulated {val['ber']:.4f} vs "
                f"theory {theory_val:.4f} ({'PASS' if passed else 'FAIL'})\n\n")
        f.write(f"{'Eb/N0(dB)':>10} {'BER_sim':>12} {'errors':>8} {'bits':>12} "
                f"{'CI_lo':>12} {'CI_hi':>12} {'BER_theory':>12}\n")
        for r in results:
            theory = float(bpsk_theory_ber(r["ebno_db"]))
            f.write(f"{r['ebno_db']:>10.1f} {r['ber']:>12.3e} {r['errors']:>8d} "
                    f"{r['bits']:>12d} {r['ci_lo']:>12.3e} {r['ci_hi']:>12.3e} "
                    f"{theory:>12.3e}\n")

    print()
    print(f"  {args.prefix}_curve.png")
    print(f"  {args.prefix}_results.txt")


if __name__ == "__main__":
    main()
