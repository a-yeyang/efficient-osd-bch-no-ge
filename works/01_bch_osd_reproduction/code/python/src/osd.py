"""Ordered Statistics Decoding (OSD) baseline for BCH codes.

Reference: Fossorier–Lin 1995 (paper's ref [5]) and the review in Sec. II-B of
the LLOSD paper.

Pipeline for order τ:
  1) Compute LLRs L_j; hard-decision r_j = (L_j < 0).
  2) Sort by |L_j| descending → permutation π.
  3) Apply π to columns of G → π(G); do Gaussian elimination so the first k
     columns form the identity submatrix; the corresponding k positions are
     the MRIPs.
  4) The (permuted) received k-bit MRIP vector is the initial message f.
  5) For each TEP e^(ω) with weight ≤ τ, form f^(ω) = f XOR e^(ω), re-encode
     ĉ^(ω) = f^(ω) · G_BCH, un-permute, and compute correlation distance
     D(r, ĉ^(ω)) = Σ_{j: r_j ≠ ĉ_j} |L_j|.
  6) Keep the candidate with minimum D. Optionally early-terminate via the ML
     criterion (eq. 14).

We also count operation counters (F2 ops, floating ops, F_2^m ops = 0 here) so
we can populate Tables I–IV directly.
"""
from __future__ import annotations

from itertools import combinations
from dataclasses import dataclass, field
from typing import Optional

import numpy as np


@dataclass
class OpCounters:
    f2: int = 0        # binary XOR / AND ops
    f2m: int = 0       # GF(2^m) multiplications / additions
    fp: int = 0        # floating-point ops (adds and comparisons)
    latency_us: float = 0.0

    def add(self, other: "OpCounters") -> "OpCounters":
        return OpCounters(
            f2=self.f2 + other.f2,
            f2m=self.f2m + other.f2m,
            fp=self.fp + other.fp,
            latency_us=self.latency_us + other.latency_us,
        )


@dataclass
class DecoderStats:
    n_frames: int = 0
    n_errors: int = 0
    ops_f2: list = field(default_factory=list)
    ops_f2m: list = field(default_factory=list)
    ops_fp: list = field(default_factory=list)
    latency_us: list = field(default_factory=list)
    n_bch_candidates: list = field(default_factory=list)
    n_teps_processed: list = field(default_factory=list)
    n_tvs_processed: list = field(default_factory=list)

    def fer(self) -> float:
        return self.n_errors / max(1, self.n_frames)


def sort_permutation_by_llr(L: np.ndarray) -> np.ndarray:
    """Return permutation π such that |L_π[0]| >= |L_π[1]| >= ..."""
    return np.argsort(-np.abs(L))


def gaussian_elim_binary(G_perm: np.ndarray):
    """Row-reduce a permuted k x n binary generator matrix so that the first k
    columns form the identity. Returns (G_sys, extra_perm) where extra_perm is
    a permutation OF COLUMNS (applied after the initial LLR permutation) which
    may be needed if the first k columns are not linearly independent. Also
    returns count of binary XOR operations performed.
    """
    G = G_perm.copy() % 2
    k, n = G.shape
    col_perm = np.arange(n)
    f2_ops = 0
    for i in range(k):
        # Find a column with a pivot at row i, starting from column i.
        piv_col = None
        for c in range(i, n):
            if G[i, c]:
                piv_col = c
                break
        if piv_col is None:
            # Need to swap a lower row up first.
            piv_row = None
            for r in range(i + 1, k):
                for c in range(i, n):
                    if G[r, c]:
                        piv_row = r
                        piv_col = c
                        break
                if piv_row is not None:
                    break
            if piv_row is None:
                raise RuntimeError("degenerate G during GE — should not happen")
            G[[i, piv_row]] = G[[piv_row, i]]
        if piv_col != i:
            # Swap columns i and piv_col; record in col_perm.
            G[:, [i, piv_col]] = G[:, [piv_col, i]]
            col_perm[[i, piv_col]] = col_perm[[piv_col, i]]
        # Eliminate this column in all other rows.
        for r in range(k):
            if r != i and G[r, i]:
                G[r] ^= G[i]
                f2_ops += n  # XOR of n bits
    return G, col_perm, f2_ops


def enumerate_teps(k: int, tau: int):
    """Yield TEPs (length-k binary vectors) with Hamming weight ≤ τ in
    increasing weight order. Weight-0 (the all-zero TEP) comes first."""
    idx = np.arange(k)
    for w in range(tau + 1):
        for support in combinations(idx.tolist(), w):
            e = np.zeros(k, dtype=np.int8)
            if w > 0:
                e[list(support)] = 1
            yield e, w


def correlation_distance(L: np.ndarray, r_hard: np.ndarray, c: np.ndarray) -> float:
    """D(r, c) = Σ_{j: r_j != c_j} |L_j|. Uses |L| directly."""
    diff = (r_hard != c)
    return float(np.abs(L[diff]).sum())


def ml_lower_bound_ok(L_sorted_abs: np.ndarray, d: int, d_omega: int, D_val: float) -> bool:
    """Eq. (14): if D(r, ĉ) ≤ sum_{j=0..d-d_ω-1} |L_ξ_j| where ξ is the
    positions where c matches r sorted by |L| ascending, then it's ML.

    In practice L_sorted_abs is |L| for the *matching* positions sorted
    ascending, and we compare D_val against the sum of the first d-d_ω-1
    of them.
    """
    K = max(0, d - d_omega - 1)
    if K == 0:
        return D_val <= 0
    return D_val <= float(L_sorted_abs[:K].sum())


def osd_decode(code, L: np.ndarray, tau: int, use_early_terminate: bool = True):
    """Full order-τ OSD. Returns (ĉ_opt, stats_dict)."""
    n = code.n
    k = code.k
    d = code.d_design
    counters = OpCounters()

    import time
    t0 = time.perf_counter()

    # Hard-decision received word from LLR sign.
    r_hard = (L < 0).astype(np.int8)
    counters.fp += n  # sign comparisons

    # Sort by |L| descending.
    perm = sort_permutation_by_llr(L)
    counters.fp += n * int(np.log2(max(n, 2)))  # rough sort cost

    # Apply permutation to G columns and to r_hard, |L|.
    Gp = code.G[:, perm]
    r_perm = r_hard[perm]
    L_perm = L[perm]

    # Gaussian elimination on Gp to expose the k x k identity in first k cols.
    G_sys, col_perm, ge_f2 = gaussian_elim_binary(Gp)
    counters.f2 += ge_f2
    # Effective column permutation seen by decoder = perm[col_perm]
    perm_eff = perm[col_perm]
    r_sorted = r_hard[perm_eff]
    L_sorted = L[perm_eff]
    absL_sorted = np.abs(L_sorted)

    # Initial message = first k bits of the sorted received word.
    f_init = r_sorted[:k].copy()

    # Precompute |L| ordering for later inequality (14) check on the matching
    # positions between r and ĉ.
    best_c_perm = None
    best_D = np.inf

    n_teps = 0
    for tep, w in enumerate_teps(k, tau):
        n_teps += 1
        f_omega = (f_init ^ tep).astype(np.int8)
        # Re-encode: c^(ω) in sorted space = f_omega @ G_sys (mod 2).
        c_sorted = (f_omega @ G_sys) % 2
        counters.f2 += k * n  # matrix-vector product

        D = correlation_distance(L_sorted, r_sorted, c_sorted)
        counters.fp += n  # sum + compare
        if D < best_D:
            best_D = D
            best_c_perm = c_sorted.copy()
            if use_early_terminate:
                # Eq. (14): compute d_ω = Hamming distance between r_sorted and c_sorted
                d_omega = int((r_sorted != c_sorted).sum())
                # ξ = positions where c matches r, sorted by |L| ascending.
                match = np.where(r_sorted == c_sorted)[0]
                if match.size > 0:
                    absL_match_sorted = np.sort(absL_sorted[match])
                    if ml_lower_bound_ok(absL_match_sorted, d, d_omega, D):
                        # Early-terminate: this is ML.
                        break

    # Un-permute the best candidate: c[perm_eff[i]] = best_c_perm[i]
    c_hat = np.zeros(n, dtype=np.int8)
    c_hat[perm_eff] = best_c_perm
    counters.latency_us = (time.perf_counter() - t0) * 1e6

    stats = {
        "counters": counters,
        "n_teps": n_teps,
        "n_bch_candidates": n_teps,  # every OSD TEP produces a BCH candidate
    }
    return c_hat, stats
