"""Numba-JIT accelerated inner loops for LLOSD/SLLOSD/HSD.

The bottleneck of the reference NumPy LLOSD implementation is the Python-level
TEP enumeration + per-TEP XOR + binary check. We rewrite the hot loop in a
Numba-compiled function that operates on the pre-built (k' × (n-k')) parity
slice of G_RS and returns the best BCH candidate (if any) plus its metadata.

The outer LLOSD wrapper (build G_RS, sort by |L|, unpermute etc.) remains
NumPy Python.

Two hot kernels:
   * order_up_to_tau_iter : enumerate all weight-<=tau TEPs, XOR G_pc rows,
     binary-check, and correlation-distance test all inline.
   * best-D tracker with ML early termination.
"""
from __future__ import annotations

import numpy as np
from numba import njit


@njit(cache=True, fastmath=False)
def llosd_inner(
    G_pc,            # (k', n-k') int64: parity slice of G_RS
    v_hat0_pc,       # (n-k',) int64: initial parity vector
    Theta,           # (k',) int64: MRIP positions in original space
    Theta_c,         # (n-k',) int64: parity positions in original space
    u0_i8,           # (k',) int8: initial message = r|_Theta
    r_hard,          # (n,) int8
    L,               # (n,) float64
    absL,            # (n,) float64
    d_design,        # int: paper's d
    tau,             # int: order
    early_terminate, # bool
):
    n = r_hard.shape[0]
    k_prime, n_minus_k_prime = G_pc.shape

    # Best candidate storage
    best_D = np.inf
    best_c = np.zeros(n, dtype=np.int8)
    have_best = False

    n_teps = 0
    n_bch = 0
    terminated_early = False

    # Weight-0 TEP: v_parity = v_hat0_pc
    for w in range(tau + 1):
        # Enumerate all combinations of `w` indices from 0..k_prime-1.
        # We do it iteratively using an index stack.
        if w == 0:
            supports_iter = [(np.zeros(0, dtype=np.int64),)]
        # else: we'll handle by nested loops below; single-flat combination
        # generator via recursion is slow in Numba, so we specialize.

    # For efficiency and to keep the JIT happy, we specialize up to tau=8.
    # We use a simple approach: enumerate via nested index arrays.
    # Below is written for tau up to 6.

    def _try(v_parity, support_size, base_c):
        pass  # placeholder; we inline below instead.

    # ---- weight 0 ----
    v_parity = v_hat0_pc.copy()
    n_teps += 1
    is_binary = True
    for jc in range(n_minus_k_prime):
        if v_parity[jc] > 1:
            is_binary = False
            break
    if is_binary:
        n_bch += 1
        c_hat = np.zeros(n, dtype=np.int8)
        for ti in range(k_prime):
            c_hat[Theta[ti]] = u0_i8[ti]
        for jc in range(n_minus_k_prime):
            c_hat[Theta_c[jc]] = np.int8(v_parity[jc])
        D = 0.0
        d_omega = 0
        for j in range(n):
            if r_hard[j] != c_hat[j]:
                D += absL[j]
                d_omega += 1
        if D < best_D:
            best_D = D
            for j in range(n):
                best_c[j] = c_hat[j]
            have_best = True
            if early_terminate:
                # ML criterion (eq. 14)
                # Collect |L| at match positions, sort ascending, sum first K.
                K = d_design - d_omega - 1
                if K <= 0:
                    if D <= 0.0:
                        terminated_early = True
                else:
                    n_match = n - d_omega
                    if n_match > 0:
                        match_L = np.empty(n_match, dtype=np.float64)
                        idx = 0
                        for j in range(n):
                            if r_hard[j] == c_hat[j]:
                                match_L[idx] = absL[j]
                                idx += 1
                        match_L.sort()
                        s = 0.0
                        for kk in range(min(K, n_match)):
                            s += match_L[kk]
                        if D <= s:
                            terminated_early = True

    if terminated_early:
        return best_c, best_D, n_teps, n_bch, have_best, terminated_early

    # ---- generic weight w in 1..tau via nested loops with early skip ----
    # We use a Python-side helper via a recursive combination generator baked
    # into per-tau loops. Numba doesn't support recursion well, so we inline
    # up to tau=6.
    idx = np.zeros(tau + 1, dtype=np.int64)
    # v_parity stack: v_stack[w] = parity vector after XOR of w rows
    v_stack = np.zeros((tau + 1, n_minus_k_prime), dtype=np.int64)
    for jc in range(n_minus_k_prime):
        v_stack[0, jc] = v_hat0_pc[jc]

    for w in range(1, tau + 1):
        # We enumerate combinations idx[0] < idx[1] < ... < idx[w-1] in
        # 0..k_prime-1 using a manual stack.
        for j in range(w):
            idx[j] = j
        # First: XOR the initial w rows into v_stack[w].
        # We rebuild v_stack levels for correctness.
        # depth = current level (0..w). At level d, v_stack[d] = XOR of
        # G_pc[idx[0..d-1]] rows over v_hat0_pc.
        # We can just recompute for each combination initially.
        depth = 0
        # Reset stack level 0.
        # (already set above)
        # Now fill up to depth = w
        while depth < w:
            # v_stack[depth+1] = v_stack[depth] XOR G_pc[idx[depth]]
            for jc in range(n_minus_k_prime):
                v_stack[depth + 1, jc] = v_stack[depth, jc] ^ G_pc[idx[depth], jc]
            depth += 1

        done = False
        while not done:
            n_teps += 1
            # Check binary
            is_binary = True
            for jc in range(n_minus_k_prime):
                if v_stack[w, jc] > 1:
                    is_binary = False
                    break
            if is_binary:
                n_bch += 1
                c_hat = np.zeros(n, dtype=np.int8)
                for ti in range(k_prime):
                    c_hat[Theta[ti]] = u0_i8[ti]
                # Flip support bits.
                for s in range(w):
                    c_hat[Theta[idx[s]]] ^= 1
                for jc in range(n_minus_k_prime):
                    c_hat[Theta_c[jc]] = np.int8(v_stack[w, jc])
                D = 0.0
                d_omega = 0
                for j in range(n):
                    if r_hard[j] != c_hat[j]:
                        D += absL[j]
                        d_omega += 1
                if D < best_D:
                    best_D = D
                    for j in range(n):
                        best_c[j] = c_hat[j]
                    have_best = True
                    if early_terminate:
                        K = d_design - d_omega - 1
                        if K <= 0:
                            if D <= 0.0:
                                terminated_early = True
                        else:
                            n_match = n - d_omega
                            if n_match > 0:
                                match_L = np.empty(n_match, dtype=np.float64)
                                p = 0
                                for j in range(n):
                                    if r_hard[j] == c_hat[j]:
                                        match_L[p] = absL[j]
                                        p += 1
                                match_L.sort()
                                s = 0.0
                                for kk in range(min(K, n_match)):
                                    s += match_L[kk]
                                if D <= s:
                                    terminated_early = True
                        if terminated_early:
                            return best_c, best_D, n_teps, n_bch, have_best, terminated_early

            # Advance combination: standard combinatorial next-combination.
            # Move deepest index; if exhausted, backtrack.
            j = w - 1
            while j >= 0 and idx[j] == k_prime - (w - j):
                j -= 1
            if j < 0:
                done = True
                break
            idx[j] += 1
            # Reset all deeper indices.
            for k in range(j + 1, w):
                idx[k] = idx[k - 1] + 1
            # Rebuild v_stack from level j upward:
            #   v_stack[j+1] = v_stack[j] XOR G_pc[idx[j]]
            for kk in range(j, w):
                for jc in range(n_minus_k_prime):
                    v_stack[kk + 1, jc] = v_stack[kk, jc] ^ G_pc[idx[kk], jc]

    return best_c, best_D, n_teps, n_bch, have_best, terminated_early
