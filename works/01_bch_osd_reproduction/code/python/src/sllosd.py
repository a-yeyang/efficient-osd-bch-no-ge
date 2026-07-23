"""SLLOSD: Segmented LLOSD (Sec. IV of the paper).

The MRIP set Θ (size k') is split into:
  - Υ : the k highest-reliability positions (top-k of MRIP)
  - Θ \ Υ : the remaining k' − k positions (still MRIP-eligible, but less
             reliable than Υ)

For each Hamming weight ρ ∈ {0..τ} of the error pattern on Υ, we allow up to
θ_ρ errors on Θ\Υ. The set E_ρ(θ_ρ) has |E_ρ(θ_ρ)| = C(k, ρ) * Σ_{ρ'=0..θ_ρ}
C(k'-k, ρ') TEPs. The full SLLOSD(θ_0, θ_1, ..., θ_τ) enumerates their union.

Example (63, 45), (θ_0, θ_1, θ_2, θ_3) = (3, 2, 1, 0):
   Σ_ρ C(k, ρ) Σ_{ρ'≤θ_ρ} C(k'-k, ρ')
 = C(45,0)*(C(12,0)+C(12,1)+C(12,2)+C(12,3))
 + C(45,1)*(C(12,0)+C(12,1)+C(12,2))
 + C(45,2)*(C(12,0)+C(12,1))
 + C(45,3)*(C(12,0))
 = 1*299 + 45*79 + 990*13 + 14190
 = 299 + 3555 + 12870 + 14190
 = 30914 (Wait, that's LLOSD(3)'s N_TEPs; the paper reports 3854 for
 SLLOSD(3,2,1,0). Let me re-check indices...)

Actually re-reading the paper: k' = |Θ| = n − 2t which for (63,45) is 57.
Υ ⊂ Θ has size k = 45. Θ\Υ has size k'−k = 57−45 = 12.

Following eq. (45) with (θ_0, θ_1, θ_2, θ_3) = (3, 2, 1, 0):
  Σ_ρ C(k, ρ) Σ_{ρ'≤θ_ρ} C(k'-k, ρ')
= C(45,0)*Σ_{ρ'≤3} C(12,ρ') + C(45,1)*Σ_{ρ'≤2} C(12,ρ') + C(45,2)*Σ_{ρ'≤1} C(12,ρ') + C(45,3)*Σ_{ρ'≤0} C(12,ρ')
= 1*(1+12+66+220) + 45*(1+12+66) + 990*(1+12) + 14190*(1)
= 299 + 3555 + 12870 + 14190
= 30914.

That is EXACTLY LLOSD(3), which is not what the paper says (3854). So my
indexing must be off. Let me re-read...

Actually the paper does say N_TEPs=30914 for LLOSD(3) on (63,45). But for
SLLOSD(3,2,1,0), N_TEPs=3854. The way this works out is θ_ρ is applied
to the complementary segment Θ\Υ ONLY when weight over Υ equals ρ (paper says
"a maximum weight of θ_ρ over Θ\Υ"), so N_TEPs should sum. My arithmetic
matches that.

Ah — the difference: the paper's "3854" applies to a DIFFERENT split! Reading
Sec IV-A more carefully: the paper's (θ_0, θ_1, θ_2, θ_3) = (3, 2, 1, 0) is
combined with k' > n/2 but they might have chosen k' larger than n−2t.

Actually re-reading Remark 1: "|Θ^c| < |Θ|" and "k' > n/2". For (63,45): if
k' = 57 (n−2t), |Θ^c| = 6, and k'-k = 12. That gives 30914 which matches
LLOSD(3). For SLLOSD to hit 3854, they must be using a *tighter* segmentation.

Looking again at eq. (45): N_TEPs = Σ_ρ C(k, ρ) Σ_{ρ'=0}^{θ_ρ} C(k'-k, ρ').
With k=45, k'-k=12, (θ_0..θ_3)=(3,2,1,0):
 ρ=0: C(45,0)*[C(12,0)+C(12,1)+C(12,2)+C(12,3)] = 1 * 299 = 299
 ρ=1: C(45,1)*[C(12,0)+C(12,1)+C(12,2)] = 45 * 79 = 3555
 ρ=2: C(45,2)*[C(12,0)+C(12,1)] = 990 * 13 = 12870
 ρ=3: C(45,3)*[C(12,0)] = 14190 * 1 = 14190
 Total = 30914.

The paper's 3854 must correspond to (θ_0..θ_3) = (?,?,?,?) that yields 3854.
Solving: 3854 = 299 + 3555 tells me only ρ=0 and ρ=1 contribute with
ρ=1's θ_1 including up to C(12,2). So actually (θ_0, θ_1) = (3, 2) and no
higher weights! I.e. τ_effective = 1 with segmentation.

But paper says "SLLOSD(3, 2, 1, 0) approximates LLOSD(3)". Let me re-read...

Hmm, in the paper: "Applying the aforementioned TEP segmentation, the LLOSD
(3) can be expressed as the SLLOSD (3, 2, 1, 0), where N_TEPs = 30914. Since
both the OSD (1) and the LLOSD (3) can approach near ML performance, it is
reasonable to expect that the SLLOSD (3, 2) can approximate the LLOSD (3),
reducing N_TEPs to 3854."

Aha! So (3, 2, 1, 0) is just re-writing LLOSD(3) in segmented form (verified:
30914). The actual reduction is SLLOSD(3, 2) — a *shorter* tuple, which means
τ = 1 in segmentation (only ρ = 0 and ρ = 1 considered).

Let's verify SLLOSD(3, 2):
 ρ=0: C(45,0)*Σ_{ρ'≤3} C(12,ρ') = 299
 ρ=1: C(45,1)*Σ_{ρ'≤2} C(12,ρ') = 45*79 = 3555
 Total = 3854. ✓ MATCHES PAPER.

So the "τ" of an SLLOSD is len(θ_tuple)−1. Great, this is clearer now.
"""
from __future__ import annotations

from itertools import combinations
import time
import numpy as np

from .osd import OpCounters, sort_permutation_by_llr, correlation_distance, ml_lower_bound_ok
from .llosd import build_rs_systematic_generator


def sllosd_decode(
    code,
    L: np.ndarray,
    theta_tuple,
    use_binary_reencoding: bool = True,
    use_early_terminate: bool = True,
):
    """Segmented LLOSD.

    theta_tuple = (θ_0, θ_1, ..., θ_τ) where τ = len(theta_tuple)-1.
    Interpretation: for each ρ ∈ 0..τ, TEPs with weight ρ on Υ (top-k of MRIP)
    combined with weight up to θ_ρ on Θ\Υ.
    """
    n = code.n
    k = code.k
    m = code.m
    d = code.d_design
    t = code.t
    gf = code.gf
    k_prime = n - 2 * t
    tau = len(theta_tuple) - 1
    counters = OpCounters()
    t0 = time.perf_counter()

    r_hard = (L < 0).astype(np.int8)
    counters.fp += n
    perm = sort_permutation_by_llr(L)
    counters.fp += n * int(np.log2(max(n, 2)))
    Theta = perm[:k_prime].copy()
    Upsilon = Theta[:k].copy()        # top-k MRIP positions
    Theta_minus_Y = Theta[k:].copy()  # remaining k' − k MRIP positions

    G_RS, Theta_c = build_rs_systematic_generator(gf, Theta, k_prime, n)
    counters.f2m += 2 * (n * n - k_prime * k_prime + k_prime)

    # Row order of G_RS matches Θ (Υ first, then Θ\Υ).
    Y_rows = np.arange(k)            # row indices for Upsilon (support of Υ)
    W_rows = np.arange(k, k_prime)   # row indices for Θ\Υ

    u0 = r_hard[Theta].astype(np.int64)
    v_hat0_pc = np.zeros(Theta_c.size, dtype=np.int64)
    active = np.where(u0 == 1)[0]
    for i in active:
        v_hat0_pc ^= G_RS[i, Theta_c]
    counters.f2m += active.size * Theta_c.size

    absL = np.abs(L)
    best_c = None
    best_D = np.inf
    n_teps = 0
    n_bch = 0

    for rho in range(tau + 1):
        theta_rho = theta_tuple[rho]
        for supp_Y in combinations(Y_rows.tolist(), rho):
            for rho_prime in range(theta_rho + 1):
                for supp_W in combinations(W_rows.tolist(), rho_prime):
                    n_teps += 1
                    support = list(supp_Y) + list(supp_W)
                    v_parity = v_hat0_pc.copy()
                    for i in support:
                        v_parity ^= G_RS[i, Theta_c]
                        if use_binary_reencoding:
                            counters.f2 += m * Theta_c.size
                        else:
                            counters.f2m += Theta_c.size
                    if np.any(v_parity > 1):
                        continue
                    n_bch += 1
                    c_hat = np.zeros(n, dtype=np.int8)
                    c_hat[Theta] = u0.astype(np.int8)
                    for i in support:
                        c_hat[Theta[i]] ^= 1
                    c_hat[Theta_c] = v_parity.astype(np.int8)
                    D = correlation_distance(L, r_hard, c_hat)
                    counters.fp += n
                    if D < best_D:
                        best_D = D
                        best_c = c_hat.copy()
                        if use_early_terminate:
                            d_omega = int((r_hard != c_hat).sum())
                            match = np.where(r_hard == c_hat)[0]
                            if match.size > 0:
                                absL_match_sorted = np.sort(absL[match])
                                if ml_lower_bound_ok(absL_match_sorted, d, d_omega, D):
                                    counters.latency_us = (time.perf_counter() - t0) * 1e6
                                    return best_c, {
                                        "counters": counters,
                                        "n_teps": n_teps,
                                        "n_bch_candidates": n_bch,
                                    }
    if best_c is None:
        best_c = r_hard.copy()
    counters.latency_us = (time.perf_counter() - t0) * 1e6
    return best_c, {
        "counters": counters,
        "n_teps": n_teps,
        "n_bch_candidates": n_bch,
    }


def sllosd_n_teps_theoretical(k: int, k_prime: int, theta_tuple) -> int:
    """Return the theoretical total N_TEPs given eq. (45)."""
    from math import comb
    total = 0
    for rho, theta_rho in enumerate(theta_tuple):
        left = comb(k, rho)
        right = sum(comb(k_prime - k, rp) for rp in range(theta_rho + 1))
        total += left * right
    return total
