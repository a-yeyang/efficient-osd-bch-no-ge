"""Maximum-Likelihood (ML) decoding via full codebook search.

For short codes (k ≤ 21 → 2M codewords) this is tractable in vectorized numpy.
For (63, 45) etc., we approximate ML by a much larger LLOSD order (τ ≥ 5) —
the paper does the same (their "ML" curves in Fig. 3/4 are effectively OSD
with τ = k, which is intractable, or use an upper bound like Union Bound
Nakagami). For (31, 21) we do full ML; for others we use an "approximate ML"
using LLOSD with a high order.
"""
from __future__ import annotations

import numpy as np


def ml_decode_full_codebook(code, L):
    """Full ML: enumerate all 2^k codewords, pick the one with minimum
    correlation distance / equivalently maximum sum of L·x."""
    k = code.k
    n = code.n
    if k > 22:
        raise ValueError(f"k={k} too large for full-codebook ML")
    # Build all codewords once as (2^k, n) matrix — cache statically.
    if not hasattr(code, "_ml_cw"):
        msgs = np.arange(2 ** k, dtype=np.int64)
        msg_bits = ((msgs[:, None] >> np.arange(k)[None, :]) & 1).astype(np.int8)
        cws = (msg_bits @ code.G) % 2
        code._ml_cw = cws
    cws = code._ml_cw
    # BPSK: bit 0 -> +1, bit 1 -> -1. Correlation = L @ (1 - 2*c).
    # Pick maximum correlation.
    corr = cws @ (-2.0 * L) + L.sum()   # equivalent to Σ L_j (1 - 2 c_j)
    best = int(np.argmax(corr))
    return cws[best], {"counters": None}


def ml_approx_by_high_order_llosd(code, L, tau=5):
    """Approximate ML by running LLOSD with a high order."""
    from .decoders import llosd_fast
    return llosd_fast(code, L, tau=tau, use_early_terminate=True)
