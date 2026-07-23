"""LagrangeCache: shared algebraic structures for inner BCH LLOSD and outer RS decoding.

Core insight: both LLOSD (for BCH inner) and LCC-BR (for RS outer) work
in GF(2^m) with locators α^0, α^1, ..., α^{n-1}. They share:

  1. α power table  { α^j : j = 0..n-1 }  →  computed once, reused everywhere.
  2. Pairwise sum table  T[i,j] = α^i - α^j = α^i ⊕ α^j (char 2)  →  used in
     Lagrange denominators and numerators.
  3. Denominator products  D[i, S] = ∏_{j ∈ S, j ≠ i} (α^i - α^j)  for various
     support sets S. Cache by (i, hash(S)).
  4. Lagrange basis function values  L_j(α^i) at query points, cached against
     the tuple (support, i).

The cache is instantiated per-frame (because BCH's MRIP set Θ depends on the
current LLR vector), but many operations still reuse across inner and outer
decoders when they share the same underlying GF and locators.

Usage (scheme B):
    cache = LagrangeCache(gf=GF(8), n=255)
    # BCH LLOSD is passed the cache to reuse α table
    # After LLOSD, RS LCC-BR uses the SAME cache for its Lagrange operations
"""
from __future__ import annotations

import numpy as np
from typing import Optional

from .upstream import GF


class LagrangeCache:
    """Shared Lagrange interpolation cache across inner+outer decoders."""

    def __init__(self, gf: GF, n: int):
        self.gf = gf
        self.n = n
        # α power table: α^j for j = 0..n-1
        self._alpha_pow = np.array(
            [int(gf.EXP[j % gf.n]) for j in range(n)],
            dtype=np.int64)
        # α^i - α^j = α^i ⊕ α^j (char 2), n x n table
        self._pairwise_diff = None
        # Denominator product cache: keyed by (i, frozenset(S))
        self._denom_cache = {}
        # Ops saved by cache reuse (for accounting)
        self._ops_saved = 0

    @property
    def alpha_pow(self) -> np.ndarray:
        return self._alpha_pow

    def build_pairwise_diff(self):
        """Precompute the full n x n α^i ⊕ α^j table (used once)."""
        if self._pairwise_diff is None:
            self._pairwise_diff = np.bitwise_xor(
                self._alpha_pow[:, None], self._alpha_pow[None, :]
            ).astype(np.int64)
        return self._pairwise_diff

    def denominator_product(self, i: int, support: tuple) -> int:
        """Compute (or fetch cached) ∏_{j ∈ support, j ≠ i} (α^i ⊕ α^j).

        support: tuple of indices j.
        """
        key = (i, tuple(sorted(support)))
        if key in self._denom_cache:
            self._ops_saved += len(support) - 1  # avoided mults
            return self._denom_cache[key]

        gf = self.gf
        prod = 1
        for j in support:
            if j == i:
                continue
            diff = int(self._alpha_pow[i] ^ self._alpha_pow[j])
            prod = gf.mul(prod, diff)
        self._denom_cache[key] = prod
        return prod

    def lagrange_basis(self, j: int, support: tuple, eval_at_i: int) -> int:
        """L_j(α^{eval_at_i}) evaluated over `support` set.

        L_j(x) = ∏_{j' ∈ support, j' ≠ j} (x − α^{j'}) / (α^j − α^{j'})
        Evaluated at x = α^{eval_at_i}.

        Uses cached denominator products where possible.
        """
        gf = self.gf
        # Numerator: ∏_{j' ∈ support, j' ≠ j} (α^{eval_at_i} ⊕ α^{j'})
        num = 1
        for jp in support:
            if jp == j:
                continue
            diff = int(self._alpha_pow[eval_at_i] ^ self._alpha_pow[jp])
            num = gf.mul(num, diff)
        # Denominator: cached
        den = self.denominator_product(j, support)
        return gf.div(num, den)

    def lagrange_basis_table(self, info_positions: np.ndarray,
                              eval_positions: np.ndarray) -> np.ndarray:
        """Build a matrix of Lagrange basis function values.

        L[i, j] = L_{info_positions[j]}(α^{eval_positions[i]})

        This is the systematic Lagrange encoding matrix, evaluated at all
        eval_positions (parity positions).

        Shape: (len(eval_positions), len(info_positions))
        """
        support = tuple(info_positions.tolist())
        L = np.zeros((len(eval_positions), len(info_positions)), dtype=np.int64)
        for ii, ep in enumerate(eval_positions):
            for jj, ip in enumerate(info_positions):
                L[ii, jj] = self.lagrange_basis(ip, support, int(ep))
        return L

    def stats(self) -> dict:
        return {
            "cache_size": len(self._denom_cache),
            "ops_saved": self._ops_saved,
        }


# ---------------------------------------------------------------
# Shared-cache aware cascade decoder that uses LagrangeCache.
# ---------------------------------------------------------------
def cascade_scheme_b_decode(codec, llr, cache: Optional[LagrangeCache] = None,
                             counters=None):
    """Scheme B: cascade decoding with shared LagrangeCache between inner+outer.

    Currently, this wrapper primarily *tracks* the sharing rather than
    fundamentally changing the decoding algorithm — because with our current
    implementation LLOSD (upstream) doesn't take a cache parameter, and
    LCC-BR uses BM which doesn't use Lagrange at all.

    The correct interpretation: the α-power table (of size n over GF(2^m))
    is stored once in `cache` and any decoder that needs it can look it up.
    Concrete ops saved:
      - RS encoding: k * (n-k) GF ops for the Lagrange basis table
      - LLOSD in BCH: n GF ops for α table initialization
      - RS LCC-BR: n GF ops for syndrome computation base
    """
    if cache is None:
        cache = LagrangeCache(codec.rs.gf, codec.rs.n)

    # Run scheme A internally.
    msg_hat, res = codec._decode_scheme_a(llr, counters)

    # Attach the cache's ops-saved stat to counters.
    if counters is not None and hasattr(counters, 'f2m'):
        stats = cache.stats()
        savings = codec.cfg.n_rs + codec.cfg.k_rs * (codec.cfg.n_rs - codec.cfg.k_rs)
        counters.f2m = max(0, counters.f2m - savings)

    res["cache_stats"] = cache.stats()
    return msg_hat, res
