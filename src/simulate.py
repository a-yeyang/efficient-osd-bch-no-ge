"""Monte Carlo simulation driver.

We keep this small and functional. Given a decoder callable, a BCH code, an
Eb/N0 sweep, and a target frame count, produce FER + latency + operation-
count statistics.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List
import time
import numpy as np

from .bch import BCHCode, bpsk_modulate, sigma_from_ebn0, llr_from_y


@dataclass
class SimulationResult:
    ebn0_db: List[float] = field(default_factory=list)
    fer: List[float] = field(default_factory=list)
    n_frames: List[int] = field(default_factory=list)
    n_errors: List[int] = field(default_factory=list)
    avg_ops_f2: List[float] = field(default_factory=list)
    avg_ops_f2m: List[float] = field(default_factory=list)
    avg_ops_fp: List[float] = field(default_factory=list)
    avg_latency_us: List[float] = field(default_factory=list)
    avg_n_bch: List[float] = field(default_factory=list)
    avg_n_teps: List[float] = field(default_factory=list)
    avg_n_tvs: List[float] = field(default_factory=list)


def run_mc(
    code: BCHCode,
    decoder,                # callable (code, L) -> (c_hat, stats_dict)
    ebn0_list,
    max_frames: int = 20_000,
    min_errors: int = 60,
    seed: int = 0,
    verbose: bool = True,
    early_stop_fer: float = 1e-6,
):
    """Run Monte Carlo simulation for one decoder over an Eb/N0 sweep.

    Continues per SNR point until either `max_frames` frames or `min_errors`
    frame errors have been observed. Uses the all-zero codeword due to
    symmetry (BPSK + AWGN + linear code).
    """
    rate = code.k / code.n
    all_zero_c = np.zeros(code.n, dtype=np.int8)
    all_zero_x = bpsk_modulate(all_zero_c)  # all +1
    res = SimulationResult()

    for ebn0 in ebn0_list:
        rng = np.random.default_rng(seed + int(ebn0 * 100))
        sigma = sigma_from_ebn0(ebn0, rate)
        n_frames = 0
        n_errors = 0
        sum_f2 = 0
        sum_f2m = 0
        sum_fp = 0
        sum_lat = 0.0
        sum_nbch = 0
        sum_nteps = 0
        sum_ntvs = 0
        t_start = time.perf_counter()
        while n_frames < max_frames:
            y = all_zero_x + sigma * rng.standard_normal(code.n)
            L = 2.0 * y / (sigma * sigma)
            c_hat, stats = decoder(code, L)
            n_frames += 1
            if not np.array_equal(c_hat, all_zero_c):
                n_errors += 1
            c = stats.get("counters")
            if c is not None:
                sum_f2 += c.f2
                sum_f2m += c.f2m
                sum_fp += c.fp
                sum_lat += c.latency_us
            sum_nbch += stats.get("n_bch_candidates", stats.get("n_bch_llosd", 0))
            sum_nteps += stats.get("n_teps", stats.get("n_teps_llosd", 0))
            sum_ntvs += stats.get("n_tvs", 0)
            if n_errors >= min_errors and n_frames >= 200:
                break

        fer = n_errors / max(1, n_frames)
        avg_f2 = sum_f2 / max(1, n_frames)
        avg_f2m = sum_f2m / max(1, n_frames)
        avg_fp = sum_fp / max(1, n_frames)
        avg_lat = sum_lat / max(1, n_frames)
        avg_nbch = sum_nbch / max(1, n_frames)
        avg_nteps = sum_nteps / max(1, n_frames)
        avg_ntvs = sum_ntvs / max(1, n_frames)
        elapsed = time.perf_counter() - t_start
        if verbose:
            print(f"  Eb/N0={ebn0:.2f} dB: FER={fer:.2e} ({n_errors}/{n_frames})  "
                  f"lat={avg_lat:.1f}us  N_BCH={avg_nbch:.2f}  N_TEPs={avg_nteps:.1f}  "
                  f"({elapsed:.1f}s)")
        res.ebn0_db.append(ebn0)
        res.fer.append(fer)
        res.n_frames.append(n_frames)
        res.n_errors.append(n_errors)
        res.avg_ops_f2.append(avg_f2)
        res.avg_ops_f2m.append(avg_f2m)
        res.avg_ops_fp.append(avg_fp)
        res.avg_latency_us.append(avg_lat)
        res.avg_n_bch.append(avg_nbch)
        res.avg_n_teps.append(avg_nteps)
        res.avg_n_tvs.append(avg_ntvs)
        if fer < early_stop_fer and n_errors >= 1:
            # Extrapolate rest; break out.
            break
    return res


def sphere_packing_bound_fer(n: int, k: int, ebn0_db_list):
    """Very rough ML approximation for BPSK-BCH: use the normal-approximation
    (NA) bound of Polyanskiy-Poor-Verdu (paper's ref [3]). We approximate ML
    performance by the union-bound-like union of pairwise error probabilities
    at the minimum distance, which is a reasonable proxy at high SNR."""
    # We use a simple Q-function based on d_min. This is a *ballpark* ML
    # curve — good for visualization, not exact.
    from math import erfc, sqrt
    def qf(x):
        return 0.5 * erfc(x / sqrt(2.0))
    # Minimum distance ≈ 2t+1. Number of nearest neighbors: A_dmin ≈
    # combinatorial. We just use A ~ k.
    # This ML curve is only for reference on the plot.
    fer = []
    for ebn0 in ebn0_db_list:
        R = k / n
        ebn0_lin = 10 ** (ebn0 / 10.0)
        # Approximate d_min as 2*ceil((n-k)/m)+1 or just n-k+1 (RS bound is
        # too loose for BCH). We just use n-k+1 as a proxy for the min dist.
        d_min = max(3, min(2 * ((n - k) // 8) + 1, n - k + 1))
        Pb = qf(sqrt(2 * R * d_min * ebn0_lin))
        # FER ≈ 1 - (1 - Pb)^n_neighbor
        fer.append(min(1.0, k * Pb))
    return fer
