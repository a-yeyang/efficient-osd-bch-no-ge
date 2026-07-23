"""Numba-JIT inner loop for SLLOSD. Same structure as llosd_inner but with
two nested weight loops (one for Upsilon = top-k of Θ, one for Θ\Υ)."""
from __future__ import annotations

import numpy as np
from numba import njit


@njit(cache=True, fastmath=False)
def sllosd_inner(
    G_pc,            # (k', n-k') int64
    v_hat0_pc,       # (n-k',) int64
    Theta,           # (k',) int64
    Theta_c,         # (n-k',) int64
    u0_i8,           # (k',) int8
    r_hard,          # (n,) int8
    L,               # (n,) float64
    absL,            # (n,) float64
    d_design,        # int
    k,               # int: BCH dimension (= size of Upsilon)
    theta_tuple,     # (tau+1,) int64: (θ_0, ..., θ_τ)
    early_terminate  # bool
):
    n = r_hard.shape[0]
    k_prime, n_minus_k_prime = G_pc.shape
    tau = theta_tuple.shape[0] - 1

    best_D = np.inf
    best_c = np.zeros(n, dtype=np.int8)
    have_best = False
    n_teps = 0
    n_bch = 0
    terminated_early = False

    # Y_rows = 0..k-1, W_rows = k..k'-1
    # For each rho in 0..tau, enumerate C(k, rho) supports of Y, then for each,
    # enumerate C(k'-k, rho') supports of W for rho' in 0..theta_tuple[rho].
    # We use manual combination enumeration via idx stacks and level-based XOR.

    # Reusable stacks: dimensions bounded by tau + max(theta) + 1
    max_theta = 0
    for i in range(tau + 1):
        if theta_tuple[i] > max_theta:
            max_theta = theta_tuple[i]
    max_depth = tau + max_theta + 1

    idx_y = np.zeros(tau + 1, dtype=np.int64)
    idx_w = np.zeros(max_theta + 1, dtype=np.int64)
    v_y = np.zeros((tau + 1, n_minus_k_prime), dtype=np.int64)
    v_w = np.zeros((max_theta + 1, n_minus_k_prime), dtype=np.int64)
    for jc in range(n_minus_k_prime):
        v_y[0, jc] = v_hat0_pc[jc]

    for rho in range(tau + 1):
        theta_rho = theta_tuple[rho]

        # Enumerate combinations of size rho from Y_rows = 0..k-1.
        def _init_idx(idx, r, cap):
            for j in range(r):
                idx[j] = j
        # Set idx_y for weight rho.
        if rho == 0:
            # Only one "empty" support.
            # v_y[0] already equals v_hat0_pc; treat as level 0.
            y_done = False
            first_y = True
        else:
            for j in range(rho):
                idx_y[j] = j
            # Fill v_y[1..rho]
            for level in range(rho):
                for jc in range(n_minus_k_prime):
                    v_y[level + 1, jc] = v_y[level, jc] ^ G_pc[idx_y[level], jc]
            first_y = True
            y_done = False

        while not y_done:
            # Current parity from Y = v_y[rho]. Now enumerate W supports.
            v_from_y_pc = v_y[rho] if rho > 0 else v_hat0_pc  # (n-k',)

            for rho_prime in range(theta_rho + 1):
                if rho_prime == 0:
                    # Empty W support; check v_from_y_pc directly.
                    n_teps += 1
                    is_binary = True
                    for jc in range(n_minus_k_prime):
                        if v_from_y_pc[jc] > 1:
                            is_binary = False
                            break
                    if is_binary:
                        n_bch += 1
                        c_hat = np.zeros(n, dtype=np.int8)
                        for ti in range(k_prime):
                            c_hat[Theta[ti]] = u0_i8[ti]
                        if rho > 0:
                            for s in range(rho):
                                c_hat[Theta[idx_y[s]]] ^= 1
                        for jc in range(n_minus_k_prime):
                            c_hat[Theta_c[jc]] = np.int8(v_from_y_pc[jc])
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
                                        sm = 0.0
                                        for kk in range(min(K, n_match)):
                                            sm += match_L[kk]
                                        if D <= sm:
                                            terminated_early = True
                                            return best_c, best_D, n_teps, n_bch, have_best, terminated_early
                    continue
                # rho_prime >= 1: W indices in [k, k'-1], size rho_prime.
                for j in range(rho_prime):
                    idx_w[j] = k + j
                # Initialize v_w[0..rho_prime] where v_w[0] = v_from_y_pc
                for jc in range(n_minus_k_prime):
                    v_w[0, jc] = v_from_y_pc[jc]
                for level in range(rho_prime):
                    for jc in range(n_minus_k_prime):
                        v_w[level + 1, jc] = v_w[level, jc] ^ G_pc[idx_w[level], jc]

                w_done = False
                while not w_done:
                    n_teps += 1
                    is_binary = True
                    for jc in range(n_minus_k_prime):
                        if v_w[rho_prime, jc] > 1:
                            is_binary = False
                            break
                    if is_binary:
                        n_bch += 1
                        c_hat = np.zeros(n, dtype=np.int8)
                        for ti in range(k_prime):
                            c_hat[Theta[ti]] = u0_i8[ti]
                        if rho > 0:
                            for s in range(rho):
                                c_hat[Theta[idx_y[s]]] ^= 1
                        for s in range(rho_prime):
                            c_hat[Theta[idx_w[s]]] ^= 1
                        for jc in range(n_minus_k_prime):
                            c_hat[Theta_c[jc]] = np.int8(v_w[rho_prime, jc])
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
                                        sm = 0.0
                                        for kk in range(min(K, n_match)):
                                            sm += match_L[kk]
                                        if D <= sm:
                                            terminated_early = True
                                            return best_c, best_D, n_teps, n_bch, have_best, terminated_early

                    # Advance idx_w (combination in [k, k'-1]).
                    j = rho_prime - 1
                    while j >= 0 and idx_w[j] == k_prime - (rho_prime - j):
                        j -= 1
                    if j < 0:
                        w_done = True
                        break
                    idx_w[j] += 1
                    for kk in range(j + 1, rho_prime):
                        idx_w[kk] = idx_w[kk - 1] + 1
                    for kk in range(j, rho_prime):
                        for jc in range(n_minus_k_prime):
                            v_w[kk + 1, jc] = v_w[kk, jc] ^ G_pc[idx_w[kk], jc]

            # Advance idx_y (combination in [0, k-1]).
            if rho == 0:
                y_done = True
                break
            j = rho - 1
            while j >= 0 and idx_y[j] == k - (rho - j):
                j -= 1
            if j < 0:
                y_done = True
                break
            idx_y[j] += 1
            for kk in range(j + 1, rho):
                idx_y[kk] = idx_y[kk - 1] + 1
            for kk in range(j, rho):
                for jc in range(n_minus_k_prime):
                    v_y[kk + 1, jc] = v_y[kk, jc] ^ G_pc[idx_y[kk], jc]

    return best_c, best_D, n_teps, n_bch, have_best, terminated_early
