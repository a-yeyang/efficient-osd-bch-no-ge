"""HSD: Hybrid Soft Decoding (LLOSD + LCC-BR).

Reference: Sec. V of the paper (Algorithm 2).

The HSD workflow:
  1) Run order-τ LLOSD. If it produces a codeword satisfying the ML criterion
     (eq. 14), terminate.
  2) Otherwise, run LCC-BR-style Chase decoding: identify the η least-reliable
     positions (LRPs), enumerate 2^η test-vectors by toggling those positions,
     and for each test-vector run a bounded-distance decoder to produce (at
     most) one BCH codeword candidate. Merge with the LLOSD list, pick the
     minimum-correlation-distance candidate.
  3) LLOSD-based skipping rule (Lemma 4): if d_H(ĉ_ω, r_ω) ≤ t for any LLOSD
     candidate ĉ_ω, skip that test-vector — Chase would have returned the same
     candidate.

Note on equivalence to LCC-BR:  Per the paper's footnote 4, LCC-BR with
interpolation multiplicity 1 has error-correction capability equal to t (the
BCH bounded-distance decoder). Hence for the *output codeword list* (and thus
the FER), running BM on each Chase test-vector is equivalent to LCC-BR.  The
complexity profile differs (F_2^m ops vs F_2), which we instrument
separately.
"""
from __future__ import annotations

from itertools import product
import time
import numpy as np

from .osd import OpCounters, sort_permutation_by_llr, correlation_distance, ml_lower_bound_ok
from .llosd import llosd_decode


def hsd_decode(
    code,
    L: np.ndarray,
    tau: int,
    eta: int,
    use_binary_reencoding: bool = True,
    use_early_terminate: bool = True,
):
    """Order-(τ, η) HSD.

    tau : order of the LLOSD primary decoder
    eta : number of least-reliable positions (LRPs) used to seed 2^η Chase
          test-vectors
    """
    n = code.n
    k = code.k
    m = code.m
    d = code.d_design
    t = code.t

    t0 = time.perf_counter()

    # Phase 1: LLOSD.
    c_llosd, s = llosd_decode(
        code, L, tau=tau,
        use_binary_reencoding=use_binary_reencoding,
        use_early_terminate=use_early_terminate,
    )
    counters = s["counters"]
    r_hard = (L < 0).astype(np.int8)
    absL = np.abs(L)

    # Best candidate so far.
    best_c = c_llosd
    best_D = correlation_distance(L, r_hard, c_llosd)
    counters.fp += n

    # Check whether LLOSD already met the ML criterion.
    d_omega = int((r_hard != c_llosd).sum())
    match = np.where(r_hard == c_llosd)[0]
    early_ml = False
    if match.size > 0 and use_early_terminate:
        absL_match_sorted = np.sort(absL[match])
        if ml_lower_bound_ok(absL_match_sorted, d, d_omega, best_D):
            early_ml = True
    n_tvs = 0
    n_tvs_skipped = 0
    if not early_ml:
        # Phase 2: LCC-BR-equivalent Chase decoding.
        perm = sort_permutation_by_llr(L)
        # η least-reliable positions.
        Psi = perm[-eta:] if eta > 0 else np.array([], dtype=np.int64)
        # Enumerate 2^η patterns.
        for bits in product([0, 1], repeat=eta):
            n_tvs += 1
            e_omega = np.zeros(n, dtype=np.int8)
            for i, b in enumerate(bits):
                if b:
                    e_omega[Psi[i]] = 1
            r_omega = (r_hard ^ e_omega).astype(np.int8)

            # Skipping rule (Lemma 4): if d_H(c_llosd, r_omega) ≤ t, skip.
            if int(np.count_nonzero(c_llosd ^ r_omega)) <= t:
                n_tvs_skipped += 1
                continue

            # Run BM on r_omega. BM's error-correction capability is t.
            dec, ok = code.bm_decode(r_omega)
            # LCC-BR complexity approximation (Lemma 12/13): C_int ≈ (n−k')^2
            # per test-vector; C_prf ≈ (n−k') GF ops.
            counters.f2m += (n - (n - 2 * t)) ** 2  # rough BR cost per TV
            counters.f2m += (n - 2 * t)             # partial root-finding cost
            if not ok:
                continue
            # d_H(dec, r_omega) ≤ t; dec is a BCH codeword.
            D = correlation_distance(L, r_hard, dec)
            counters.fp += n
            if D < best_D:
                best_D = D
                best_c = dec
                # Check ML.
                d_omega = int((r_hard != dec).sum())
                match = np.where(r_hard == dec)[0]
                if match.size > 0 and use_early_terminate:
                    absL_match_sorted = np.sort(absL[match])
                    if ml_lower_bound_ok(absL_match_sorted, d, d_omega, D):
                        break

    counters.latency_us = (time.perf_counter() - t0) * 1e6
    return best_c, {
        "counters": counters,
        "n_teps_llosd": s["n_teps"],
        "n_bch_llosd": s["n_bch_candidates"],
        "n_tvs": n_tvs,
        "n_tvs_skipped": n_tvs_skipped,
    }
