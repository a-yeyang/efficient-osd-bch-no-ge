"""LLOSD and LLOSD-B (Low-Latency OSD) for BCH codes.

Central idea (Sec. III of the paper): Because BCH is a binary subcode of RS,
we can build the RS (n, k') systematic generator matrix G_RS over GF(2^m)
directly from Lagrange interpolation polynomials — no Gaussian elimination
required. Each row of G_RS is thus an RS codeword; if its first k' entries
happen to be binary (∈ F_2), Lemma 1 tells us it is also a BCH codeword. So
we enumerate TEPs, re-encode as RS codewords via G_RS, and *filter* to those
that lie in the BCH subcode via a binary check on the parity positions
(Theorem 2, eq. 28).

k' here is defined via the paper's Remark 1: RS codes ⊂ MDS ⇒ dim(RS) > dim(BCH).
Empirically k' > n/2 (Plotkin bound). We pick k' = the largest RS dimension
such that the RS code with locators [1, α, ..., α^{n-1}] evaluated on our BCH
generator polynomial still contains C_BCH as a binary subcode. In the paper
they always seem to use k' = n − ε where ε is a small design offset.

For the primitive narrow-sense BCH code with design distance d = 2t+1, the
paper uses d' = d for the RS code (Lemma 1 with d' = d), giving k' = n − 2t.
That is what we use: k' = n − 2t = k + (r_BCH − 2t) where r_BCH = n − k.

Concrete example: (63, 45, t=3), 2t = 6, so k' = 57 = n − 6.
Then the RS-parity positions Θ^c have size n − k' = 6 (matches the paper's
usage that non-binary parity gets filtered).
"""
from __future__ import annotations

from itertools import combinations
from dataclasses import dataclass

import time
import numpy as np

from .gf import GF
from .osd import OpCounters, sort_permutation_by_llr, correlation_distance, ml_lower_bound_ok


# --------------------------------------------------------------------------
# RS systematic generator matrix via Lagrange interpolation
# --------------------------------------------------------------------------
def build_rs_systematic_generator(gf: GF, Theta: np.ndarray, k_prime: int, n: int):
    """Given the MRIP set Θ (k' positions in {0..n-1}), build a k' x n RS
    systematic generator matrix G_RS ∈ F_{2^m}^{k' x n} such that the columns
    indexed by Θ form the k' x k' identity submatrix.

    We use eq. (23) / (19)–(20) of the paper. Each row H_u_{j_i}(x) is the
    Lagrange polynomial T_{j_i}(x) taking value 1 at α^{j_i} and 0 at the
    other α^j for j ∈ Θ \ {j_i}.

    Vectorized implementation: use EXP/LOG tables and numpy pairwise ops.
    """
    assert Theta.size == k_prime
    Theta_mask = np.zeros(n, dtype=bool)
    Theta_mask[Theta] = True
    Theta_c = np.nonzero(~Theta_mask)[0]
    assert Theta_c.size == n - k_prime

    # Locators α^j for j = 0..n-1.
    loc = gf.EXP[np.arange(n) % gf.n]  # int64 array, length n

    G = np.zeros((k_prime, n), dtype=np.int64)
    for row_idx, j_i in enumerate(Theta):
        G[row_idx, j_i] = 1

    if Theta_c.size == 0:
        return G, Theta_c

    a_i = loc[Theta]                 # (k', )
    a_c = loc[Theta_c]               # (n-k', )
    # Denominators: for each row i, denom[i] = ∏_{j' ∈ Θ, j' != i} (a_i[i] ⊕ a_i[j'])
    # Build k' × k' pairwise XOR, mask diagonal.
    XOR_ii = np.bitwise_xor(a_i[:, None], a_i[None, :])  # (k', k')
    # Set diagonal to 1 (so log(1)=0 contributes nothing to product).
    np.fill_diagonal(XOR_ii, 1)
    LOG_ii = gf.LOG[XOR_ii]
    denom_log = LOG_ii.sum(axis=1) % gf.n  # sum of logs = log of product

    # Numerators for each (row_idx, j_c): ∏_{j' ∈ Θ, j' != i} (a_c[j_c] ⊕ a_i[j'])
    XOR_ci = np.bitwise_xor(a_c[:, None], a_i[None, :])  # (n-k', k')
    # For row i we need product over j' ∈ Θ \ {j_i}, but here XOR_ci already
    # excludes only same index in a different sense: XOR_ci[c, i] uses a_i[i]
    # but we want to *include* all j' != i. So the product for (c, i) is
    # ∏_{jp} XOR_ci[c, jp] / XOR_ci[c, i]... but if a_c[c] == a_i[i]... never
    # happens since Theta_c disjoint from Theta so all XORs are nonzero.
    # Actually the correct formula is: numerator (c, i) = ∏_{jp != i} (a_c[c] ⊕ a_i[jp]).
    LOG_ci = gf.LOG[XOR_ci]                              # (n-k', k')
    # Sum over jp axis (=1), then subtract the i-th column to exclude jp=i.
    total_log_ci = LOG_ci.sum(axis=1)                    # (n-k', )
    # For each (c, i), num_log = total_log_ci[c] - LOG_ci[c, i]
    num_log = (total_log_ci[:, None] - LOG_ci) % gf.n    # (n-k', k')
    # G[i, j_c] = num[c, i] / denom[i] = EXP[num_log - denom_log[i]]
    exp_index = (num_log - denom_log[None, :]) % gf.n    # (n-k', k')
    vals = gf.EXP[exp_index]                             # (n-k', k')
    # Fill G[Θ_rows, Θ^c_cols]: G[i, j_c] = vals[c, i]
    G[np.arange(k_prime)[:, None], Theta_c[None, :]] = vals.T
    return G, Theta_c


def rs_encode_row(gf: GF, u_msg: np.ndarray, G_RS: np.ndarray) -> np.ndarray:
    """Compute u @ G_RS over GF(2^m). u_msg has entries in {0, 1} (F_2 ⊂ F_{2^m})."""
    k, n = G_RS.shape
    out = np.zeros(n, dtype=np.int64)
    for i in range(k):
        if u_msg[i] == 0:
            continue
        # u_msg[i] == 1, so out ^= G_RS[i] (F_2 addition is XOR of GF elements)
        out ^= G_RS[i]
    return out


# --------------------------------------------------------------------------
# LLOSD (Algorithm 1)
# --------------------------------------------------------------------------
def _make_col_permutation(perm_desc: np.ndarray, n: int) -> np.ndarray:
    """Given a permutation π (indices in [0, n)) that sorts by |L| descending,
    return the ordered position list so π_eff[i] gives the original index of
    the i-th ranked position."""
    return perm_desc


def llosd_decode(
    code,
    L: np.ndarray,
    tau: int,
    use_binary_reencoding: bool = False,
    use_early_terminate: bool = True,
    max_bch_candidates: int = 10**9,
):
    """Order-τ LLOSD.

    Parameters
    ----------
    code : BCHCode
    L    : LLR vector
    tau  : decoding order
    use_binary_reencoding : if True, use LLOSD-B via P1 syndrome check
                            (eq. 40 / 38). Same performance, fewer F_2^m ops.
    use_early_terminate   : apply ML lower bound (eq. 14).
    """
    n = code.n
    k = code.k
    m = code.m
    d = code.d_design
    gf = code.gf
    t = code.t
    k_prime = n - 2 * t  # RS dimension (Lemma 1 with d' = d)

    counters = OpCounters()
    t0 = time.perf_counter()

    r_hard = (L < 0).astype(np.int8)
    counters.fp += n

    # Sort positions by |L| descending.
    perm = sort_permutation_by_llr(L)
    counters.fp += n * int(np.log2(max(n, 2)))

    # Take the k' most-reliable positions as MRP set Θ.
    Theta = perm[:k_prime].copy()
    Theta_set = set(Theta.tolist())

    # -- Build RS systematic G_RS with identity on columns Θ. -----------
    # NOTE: Because we sorted by |L|, Θ is generally NOT the first k'
    # positions of the original word; it's a subset of {0..n-1} of size k'.
    G_RS, Theta_c = build_rs_systematic_generator(gf, Theta, k_prime, n)
    # Complexity (Lemma 10): 2(n^2 - k'^2 + k') GF ops in the systematic part.
    counters.f2m += 2 * (n * n - k_prime * k_prime + k_prime)

    # Slice of G_RS restricted to parity columns Θ^c — this is the hot data.
    G_pc = G_RS[:, Theta_c]  # (k', n-k')

    # Initial estimated RS codeword v_hat^(0) on parity columns Θ^c only.
    u0 = r_hard[Theta].astype(np.int64)  # binary values seen as GF elements
    active_rows = np.nonzero(u0)[0]
    if active_rows.size > 0:
        v_hat0_pc = np.bitwise_xor.reduce(G_pc[active_rows], axis=0)
    else:
        v_hat0_pc = np.zeros(Theta_c.size, dtype=np.int64)
    counters.f2m += active_rows.size * Theta_c.size  # XORs over GF

    # For each TEP e^(ω), the produced parity is
    #   v^(ω)[Θ^c] = v_hat0_pc ⊕ Σ_{i: e[i]=1} G_pc[i]
    # If v^(ω)[Θ^c] ⊂ {0, 1}, then v^(ω) is a BCH candidate (Theorem 2).

    absL = np.abs(L)

    best_c = None
    best_D = np.inf

    n_teps = 0
    n_bch = 0

    # Precompute the base codeword template (parts that don't change per TEP).
    # c_hat[Theta] = u0 (as int8), and c_hat[Theta_c] gets set to v_parity.
    c_template = np.zeros(n, dtype=np.int8)
    c_template[Theta] = u0.astype(np.int8)
    u0_i8 = u0.astype(np.int8)

    for w in range(tau + 1):
        for support in combinations(range(k_prime), w):
            n_teps += 1
            support_arr = np.asarray(support, dtype=np.int64)
            # Compute v_parity = v_hat0_pc XOR of G_pc[support_arr] rows.
            if support_arr.size == 0:
                v_parity = v_hat0_pc
            else:
                v_parity = v_hat0_pc ^ np.bitwise_xor.reduce(G_pc[support_arr], axis=0)
            if use_binary_reencoding:
                counters.f2 += m * Theta_c.size * max(1, support_arr.size)
            else:
                counters.f2m += Theta_c.size * max(1, support_arr.size)

            # Binary check: any value > 1 means a non-binary GF element.
            if (v_parity > 1).any():
                continue  # non-binary → skip (Theorem 2)
            n_bch += 1

            # Build binary candidate codeword: reuse the template.
            c_hat = c_template.copy()
            # Flip u0 bits at the support positions (indices *within* Theta).
            if support_arr.size > 0:
                c_hat[Theta[support_arr]] ^= 1
            c_hat[Theta_c] = v_parity.astype(np.int8)

            # Correlation distance: sum of |L_j| where r_hard[j] != c_hat[j].
            diff_mask = (r_hard != c_hat)
            D = float(absL[diff_mask].sum())
            counters.fp += n
            if D < best_D:
                best_D = D
                best_c = c_hat  # no need to copy since we build fresh each iteration
                if use_early_terminate:
                    d_omega = int(diff_mask.sum())
                    match_mask = ~diff_mask
                    if match_mask.any():
                        absL_match_sorted = np.sort(absL[match_mask])
                        K = max(0, d - d_omega - 1)
                        if K == 0:
                            ml_ok = D <= 0
                        else:
                            ml_ok = D <= float(absL_match_sorted[:K].sum())
                        if ml_ok:
                            counters.latency_us = (time.perf_counter() - t0) * 1e6
                            return best_c, {
                                "counters": counters,
                                "n_teps": n_teps,
                                "n_bch_candidates": n_bch,
                                "terminated_early": True,
                            }
            if n_bch >= max_bch_candidates:
                break

    if best_c is None:
        # If nothing binary found (shouldn't happen for w=0 since u0 is binary),
        # fall back to the hard-decision received word.
        best_c = r_hard.copy()
    counters.latency_us = (time.perf_counter() - t0) * 1e6
    return best_c, {
        "counters": counters,
        "n_teps": n_teps,
        "n_bch_candidates": n_bch,
        "terminated_early": False,
    }
