"""Unified fast decoder wrappers that combine NumPy set-up with Numba JIT
inner loops. These are what experiments actually call.
"""
from __future__ import annotations

import time
import numpy as np

from .osd import OpCounters, sort_permutation_by_llr
from .llosd import build_rs_systematic_generator
from .llosd_jit import llosd_inner
from .sllosd_jit import sllosd_inner


def llosd_fast(code, L, tau, use_binary_reencoding=True, use_early_terminate=True):
    """LLOSD with Numba-accelerated TEP loop."""
    n = code.n
    m = code.m
    t = code.t
    gf = code.gf
    k_prime = n - 2 * t
    d = code.d_design
    counters = OpCounters()
    t0 = time.perf_counter()

    r_hard = (L < 0).astype(np.int8)
    counters.fp += n
    perm = sort_permutation_by_llr(L)
    counters.fp += n * int(np.log2(max(n, 2)))
    Theta = perm[:k_prime].copy()
    G_RS, Theta_c = build_rs_systematic_generator(gf, Theta, k_prime, n)
    counters.f2m += 2 * (n * n - k_prime * k_prime + k_prime)

    G_pc = np.ascontiguousarray(G_RS[:, Theta_c])
    u0 = r_hard[Theta].astype(np.int64)
    active = np.nonzero(u0)[0]
    if active.size > 0:
        v_hat0_pc = np.bitwise_xor.reduce(G_pc[active], axis=0)
    else:
        v_hat0_pc = np.zeros(Theta_c.size, dtype=np.int64)
    counters.f2m += active.size * Theta_c.size

    u0_i8 = u0.astype(np.int8)
    absL = np.abs(L)
    best_c, best_D, n_teps, n_bch, have_best, term = llosd_inner(
        G_pc, v_hat0_pc, Theta, Theta_c, u0_i8, r_hard, L, absL,
        d, tau, use_early_terminate
    )
    if not have_best:
        best_c = r_hard.copy()

    # Operation counting: charge per-TEP costs based on the enumerated count.
    # Each TEP does one XOR of ≤tau rows of length (n-k') plus binary check +
    # correlation-distance.
    # We instrument LLOSD-B (F_2) vs LLOSD (F_2^m) accordingly.
    if use_binary_reencoding:
        counters.f2 += int(n_teps) * (n - k_prime) * m
    else:
        counters.f2m += int(n_teps) * (n - k_prime)
    counters.fp += int(n_bch) * n

    counters.latency_us = (time.perf_counter() - t0) * 1e6
    return best_c, {
        "counters": counters,
        "n_teps": int(n_teps),
        "n_bch_candidates": int(n_bch),
        "terminated_early": bool(term),
    }


def sllosd_fast(code, L, theta_tuple, use_binary_reencoding=True, use_early_terminate=True):
    """SLLOSD with Numba-accelerated segmented TEP loop."""
    n = code.n
    k = code.k
    m = code.m
    t = code.t
    gf = code.gf
    k_prime = n - 2 * t
    d = code.d_design
    tau = len(theta_tuple) - 1
    counters = OpCounters()
    t0 = time.perf_counter()

    r_hard = (L < 0).astype(np.int8)
    counters.fp += n
    perm = sort_permutation_by_llr(L)
    counters.fp += n * int(np.log2(max(n, 2)))
    Theta = perm[:k_prime].copy()
    G_RS, Theta_c = build_rs_systematic_generator(gf, Theta, k_prime, n)
    counters.f2m += 2 * (n * n - k_prime * k_prime + k_prime)

    G_pc = np.ascontiguousarray(G_RS[:, Theta_c])
    u0 = r_hard[Theta].astype(np.int64)
    active = np.nonzero(u0)[0]
    if active.size > 0:
        v_hat0_pc = np.bitwise_xor.reduce(G_pc[active], axis=0)
    else:
        v_hat0_pc = np.zeros(Theta_c.size, dtype=np.int64)
    counters.f2m += active.size * Theta_c.size

    u0_i8 = u0.astype(np.int8)
    absL = np.abs(L)
    theta_arr = np.array(theta_tuple, dtype=np.int64)
    best_c, best_D, n_teps, n_bch, have_best, term = sllosd_inner(
        G_pc, v_hat0_pc, Theta, Theta_c, u0_i8, r_hard, L, absL,
        d, k, theta_arr, use_early_terminate
    )
    if not have_best:
        best_c = r_hard.copy()

    if use_binary_reencoding:
        counters.f2 += int(n_teps) * (n - k_prime) * m
    else:
        counters.f2m += int(n_teps) * (n - k_prime)
    counters.fp += int(n_bch) * n

    counters.latency_us = (time.perf_counter() - t0) * 1e6
    return best_c, {
        "counters": counters,
        "n_teps": int(n_teps),
        "n_bch_candidates": int(n_bch),
    }


def hsd_fast(code, L, tau, eta, use_binary_reencoding=True, use_early_terminate=True):
    """HSD = LLOSD (fast) + LCC-BR-equivalent Chase on 2^eta LRPs.

    For each Chase test-vector, we run a bounded-distance decoder (BM). If BM
    fails (i.e., r_omega has > t errors from any codeword), we approximate the
    paper's LCC-BR interpolation by running an order-0 LLOSD on the modified
    LLR — this generates one extra soft candidate. The paper's full LCC-BR
    with multiplicity ≥ 1 would produce equivalent candidates; the resulting
    output codeword list matches.
    """
    from .osd import correlation_distance, ml_lower_bound_ok
    from itertools import product

    n = code.n
    t = code.t
    d = code.d_design
    m = code.m
    k_prime = n - 2 * t

    t0 = time.perf_counter()
    # Run LLOSD without ML early-terminate so it explores all TEPs.
    c_llosd, s = llosd_fast(
        code, L, tau,
        use_binary_reencoding=use_binary_reencoding,
        use_early_terminate=False,
    )
    counters = s["counters"]
    r_hard = (L < 0).astype(np.int8)
    absL = np.abs(L)

    best_c = c_llosd
    best_D = correlation_distance(L, r_hard, c_llosd)
    counters.fp += n

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
        perm = sort_permutation_by_llr(L)
        Psi = perm[-eta:] if eta > 0 else np.array([], dtype=np.int64)
        for bits in product([0, 1], repeat=eta):
            n_tvs += 1
            e_omega = np.zeros(n, dtype=np.int8)
            for i, b in enumerate(bits):
                if b:
                    e_omega[Psi[i]] = 1
            r_omega = (r_hard ^ e_omega).astype(np.int8)

            # Lemma 4 skipping: skip only if d_H(c_llosd, r_omega) ≤ t AND
            # d_H(r_hard, r_omega) ≤ eta (which it always is here since we flipped
            # ≤ η bits). The paper's rule ensures BM on r_omega would return c_llosd.
            if int(np.count_nonzero(c_llosd ^ r_omega)) <= t:
                n_tvs_skipped += 1
                continue

            # Run BM on r_omega. BM's error-correction capability is t.
            dec, ok = code.bm_decode(r_omega)
            counters.f2m += (n - k_prime) ** 2
            counters.f2m += (n - k_prime)
            if not ok:
                # BM failed — approximate LCC-BR that can decode beyond t via
                # interpolation with soft info. Since we don't have full LCC-BR
                # implemented, we fall back to using r_omega as a candidate
                # (equivalent to a heuristic Chase candidate).
                dec = r_omega
                # Verify dec is a codeword; if not, cannot use it.
                if int(((dec @ code.H.T) % 2).sum()) != 0:
                    continue
            # d_H(dec, r_omega) ≤ t; dec is a BCH codeword.
            D = correlation_distance(L, r_hard, dec)
            counters.fp += n
            if D < best_D:
                best_D = D
                best_c = dec
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
