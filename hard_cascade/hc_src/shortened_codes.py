"""Shortened RS and binary BCH codes for the v3 cascade update.

Motivation (advisor's third task, 2026-07-27):
  Config 1 (NEW):  RS(544,514, t=15) over GF(2^10)  +  BCH(144,136, t=1) over GF(2^8)
  Config 2 (KEEP): RS(255,239, t=8)  over GF(2^8)   +  BCH(255,239, t=2) over GF(2^8)

Neither RS(544,514) nor BCH(144,136) is a natural (n = 2^m - 1) code, so both
are *shortened* codes.  We wrap the tested full-length decoders:

  * ShortenedRSCode(m, n_short, k_short)  wraps RSCode(m, k_full)
  * BinaryBCH(m, t, n_short, k_short)     generalises BCHt2Code to t in {1, 2}

Shortening principle
--------------------
A systematic codeword has layout  c = [ parity (n-k) | message (k) ].
To shorten by  s = k_full - k_short  symbols we force the first ``s`` message
symbols to zero.  Those zeros are never transmitted and are *known* to the
decoder, so they carry no channel noise:

  encode:  full_msg = [0]*s ++ real_msg           (length k_full)
           c_full    = [ parity(n-k) | 0*s | real_msg ]
           tx        = [ parity(n-k) | real_msg ]  (length n_short = n_full - s)

  decode:  r_full = [ rx_parity | 0*s | rx_msg ]   (re-insert known zeros)
           run full-length decoder  ->  c_dec_full
           real_msg_hat = c_dec_full[(n-k)+s : ]

Because the padding zeros are noiseless, the full-length BM / Chien / Forney
math is unchanged; error positions returned in full coordinates are correct.
"""
from __future__ import annotations

import numpy as np

from .upstream import GF, OpCounters, RSCode


# =====================================================================
# Shortened Reed-Solomon
# =====================================================================
class ShortenedRSCode:
    """Shortened narrow-sense RS(n_short, k_short) over GF(2^m).

    Built from the full RS(N=2^m-1, K_full) with the same parity count
    (n_short - k_short == N - K_full == 2t).
    """

    def __init__(self, m: int, n_short: int, k_short: int):
        self.m = m
        self.n = n_short
        self.k = k_short
        self.n_parity = n_short - k_short
        assert self.n_parity % 2 == 0, "RS parity must be even (=2t)"
        self.t = self.n_parity // 2

        self.N_full = (1 << m) - 1
        self.K_full = self.N_full - self.n_parity
        self.shorten = self.K_full - self.k          # s
        assert self.shorten >= 0, \
            f"cannot shorten: k_short={k_short} > K_full={self.K_full}"
        assert self.N_full - self.shorten == self.n

        self.full = RSCode(m=m, k=self.K_full)
        self.gf = self.full.gf

    def describe(self) -> str:
        tag = "" if self.shorten == 0 else f" [shortened by {self.shorten}]"
        return f"RS({self.n},{self.k}, t={self.t}) /GF(2^{self.m}){tag}"

    # ------------------------------------------------------------------
    def encode_systematic(self, msg: np.ndarray) -> np.ndarray:
        assert msg.size == self.k
        full_msg = np.concatenate(
            [np.zeros(self.shorten, dtype=np.int64), msg.astype(np.int64)])
        c_full = self.full.encode_systematic(full_msg)
        # c_full = [ parity(n_parity) | 0*s | real_msg ]
        parity = c_full[:self.n_parity]
        real = c_full[self.n_parity + self.shorten:]
        return np.concatenate([parity, real]).astype(np.int64)

    def extract_message(self, c_short: np.ndarray) -> np.ndarray:
        return c_short[self.n_parity:].copy()

    def _reinsert(self, r_short: np.ndarray) -> np.ndarray:
        return np.concatenate([
            r_short[:self.n_parity].astype(np.int64),
            np.zeros(self.shorten, dtype=np.int64),
            r_short[self.n_parity:].astype(np.int64),
        ])

    def bm_decode(self, r_short: np.ndarray, counters: OpCounters = None):
        r_full = self._reinsert(r_short)
        c_dec_full, ok = self.full.bm_decode(r_full, counters)
        # Drop the known padding zeros → shortened codeword coordinates
        c_short = np.concatenate([
            c_dec_full[:self.n_parity],
            c_dec_full[self.n_parity + self.shorten:],
        ]).astype(np.int64)
        return c_short, ok


# =====================================================================
# Binary BCH (t = 1 or 2), optionally shortened
# =====================================================================
class BinaryBCH:
    """Binary primitive narrow-sense BCH code, t in {1, 2}, optional shortening.

    Full length N = 2^m - 1, K_full = N - m*t.  A shortened instance
    (n_short, k_short) keeps the same parity count (m*t).

    Two hard-decision decoders:
      * decode_conventional : BM(+Chien for t=2), classic path
      * decode_direct       : Lagendijk 2026 direct root finding
                              (for t=1 both paths reduce to p = log(S1))
    """

    def __init__(self, m: int, t: int, n_short: int = None, k_short: int = None):
        assert t in (1, 2), "BinaryBCH supports t=1 and t=2"
        self.m = m
        self.t = t
        self.gf = GF(m)
        self.N_full = self.gf.n              # 2^m - 1
        self.d_design = 2 * t + 1

        # Generator polynomial g = prod of minimal polys of odd powers 1,3,..,2t-1
        seen = set()
        g_poly = [1]
        for i in range(1, 2 * t + 1):
            coset = self._cyclotomic_coset(i)
            rep = min(coset)
            if rep in seen:
                continue
            seen.add(rep)
            m_poly = [1]
            for s in coset:
                root = int(self.gf.EXP[s])
                m_poly = self.gf.poly_mul(m_poly, [root, 1])
            m_poly_bin = [int(c) & 1 for c in m_poly]
            g_poly = self.gf.poly_mul(g_poly, m_poly_bin)
            g_poly = [int(c) & 1 for c in g_poly]
        self.g_poly = g_poly
        self.n_parity = len(g_poly) - 1
        self.K_full = self.N_full - self.n_parity
        assert self.n_parity == m * t, \
            f"deg(g)={self.n_parity} != m*t={m*t}"

        # Shortening
        if n_short is None:
            n_short = self.N_full
        if k_short is None:
            k_short = self.K_full
        self.n = n_short
        self.k = k_short
        self.shorten = self.K_full - self.k
        assert self.shorten >= 0 and self.N_full - self.shorten == self.n, \
            f"bad shortening: n={n_short}, k={k_short}, K_full={self.K_full}"
        assert self.n - self.k == self.n_parity

        self.G = self._build_G()               # full-length k_full x N_full
        if t == 2:
            self.lut_A, self.lut_A_valid = self._build_lut_A()

    # ------------------------------------------------------------------
    def _cyclotomic_coset(self, i):
        coset = set()
        j = i
        while j not in coset:
            coset.add(j)
            j = (j * 2) % self.N_full
        return coset

    def _build_G(self):
        N, K = self.N_full, self.K_full
        n_minus_k = self.n_parity
        G = np.zeros((K, N), dtype=np.int8)
        gf = self.gf
        g = list(self.g_poly)
        for i in range(K):
            dividend = [0] * (n_minus_k + i) + [1] + [0] * (K - i - 1)
            dividend = dividend[:N]
            _q, rem = gf.poly_divmod(dividend, g)
            parity = list(rem) + [0] * (n_minus_k - len(rem))
            G[i, :n_minus_k] = parity[:n_minus_k]
            G[i, n_minus_k + i] = 1
        return G

    def describe(self) -> str:
        tag = "" if self.shorten == 0 else f" [shortened by {self.shorten}]"
        return f"BCH({self.n},{self.k}, t={self.t}) /GF(2^{self.m}){tag}"

    # ------------------------------------------------------------------
    # Encode / extract (shortened interface)
    # ------------------------------------------------------------------
    def encode(self, msg: np.ndarray) -> np.ndarray:
        assert msg.size == self.k
        full_msg = np.concatenate(
            [np.zeros(self.shorten, dtype=np.int64), msg.astype(np.int64)])
        c_full = (full_msg @ self.G) % 2      # [ parity | 0*s | real_msg ]
        parity = c_full[:self.n_parity]
        real = c_full[self.n_parity + self.shorten:]
        return np.concatenate([parity, real]).astype(np.int8)

    def extract_message(self, cw_short: np.ndarray) -> np.ndarray:
        return cw_short[self.n_parity:]

    def _reinsert(self, r_short: np.ndarray) -> np.ndarray:
        return np.concatenate([
            r_short[:self.n_parity].astype(np.int8),
            np.zeros(self.shorten, dtype=np.int8),
            r_short[self.n_parity:].astype(np.int8),
        ])

    def _drop(self, cw_full: np.ndarray) -> np.ndarray:
        return np.concatenate([
            cw_full[:self.n_parity],
            cw_full[self.n_parity + self.shorten:],
        ]).astype(np.int8)

    # ------------------------------------------------------------------
    # Syndromes (full-length coordinates)
    # ------------------------------------------------------------------
    def compute_syndromes(self, r_full: np.ndarray, counters: OpCounters = None):
        if counters is None:
            counters = OpCounters()
        gf = self.gf
        nz = np.nonzero(r_full)[0]
        if nz.size == 0:
            return 0, 0
        s1 = int(np.bitwise_xor.reduce(gf.EXP[nz % gf.n]))
        s3 = 0
        if self.t >= 2:
            s3 = int(np.bitwise_xor.reduce(gf.EXP[(3 * nz) % gf.n]))
        counters.f2m += (2 if self.t >= 2 else 1) * self.n
        return s1, s3

    # ------------------------------------------------------------------
    # Decoder 1: Conventional
    # ------------------------------------------------------------------
    def decode_conventional(self, r_short: np.ndarray, counters: OpCounters = None):
        gf = self.gf
        n = self.N_full
        if counters is None:
            counters = OpCounters()
        r_full = self._reinsert(r_short)
        S1, S3 = self.compute_syndromes(r_full, counters)

        if self.t == 1:
            if S1 == 0:
                return self._drop(r_full), True
            p = int(gf.LOG[S1]) % n
            c = r_full.copy()
            c[p] ^= 1
            return self._drop(c), True

        # t == 2
        S1_cubed = gf.mul(gf.mul(S1, S1), S1)
        counters.f2m += 2
        if S1 == 0 and S3 == 0:
            return self._drop(r_full), True
        if S1 == 0 and S3 != 0:
            return self._drop(r_full), False
        if S1_cubed == S3:
            p = int(gf.LOG[S1]) % n
            c = r_full.copy()
            c[p] ^= 1
            return self._drop(c), True
        c0 = gf.div(gf.add(S1_cubed, S3), S1)
        counters.f2m += 2
        err_positions = []
        for i in range(n):
            X = int(gf.EXP[i])
            X_sq = int(gf.EXP[(2 * i) % n])
            term = X_sq ^ gf.mul(S1, X) ^ c0
            counters.f2m += 2
            if term == 0:
                err_positions.append(i)
                if len(err_positions) == 2:
                    break
        if len(err_positions) != 2:
            return self._drop(r_full), False
        c = r_full.copy()
        for p in err_positions:
            c[p] ^= 1
        return self._drop(c), True

    # ------------------------------------------------------------------
    # Decoder 2: Direct root finding (Lagendijk §III-A)
    # ------------------------------------------------------------------
    def _build_lut_A(self):
        gf = self.gf
        n_field = 1 << self.m
        lut = np.zeros((n_field, 2), dtype=np.int64)
        lut_valid = np.zeros(n_field, dtype=bool)
        for X in range(n_field):
            X2 = gf.mul(X, X)
            k = X2 ^ X
            if not lut_valid[k]:
                lut[k, 0] = X
                lut[k, 1] = X ^ 1
                lut_valid[k] = True
        return lut, lut_valid

    def decode_direct(self, r_short: np.ndarray, counters: OpCounters = None):
        gf = self.gf
        n = self.N_full
        if counters is None:
            counters = OpCounters()
        r_full = self._reinsert(r_short)
        S1, S3 = self.compute_syndromes(r_full, counters)

        if self.t == 1:
            # Direct for t=1: no quadratic, no LUT — closed form p = log(S1)
            if S1 == 0:
                return self._drop(r_full), True
            p = int(gf.LOG[S1]) % n
            c = r_full.copy()
            c[p] ^= 1
            return self._drop(c), True

        # t == 2
        S1_cubed = gf.mul(gf.mul(S1, S1), S1)
        counters.f2m += 2
        if S1 == 0 and S3 == 0:
            return self._drop(r_full), True
        if S1 == 0 and S3 != 0:
            return self._drop(r_full), False
        if S1_cubed == S3:
            p = int(gf.LOG[S1]) % n
            c = r_full.copy()
            c[p] ^= 1
            return self._drop(c), True
        numerator = gf.add(S1_cubed, S3)
        k_lut = gf.div(numerator, S1_cubed)
        counters.f2m += 1
        if not self.lut_A_valid[k_lut]:
            return self._drop(r_full), False
        Y1 = int(self.lut_A[k_lut, 0])
        Y2 = int(self.lut_A[k_lut, 1])
        X1 = gf.mul(S1, Y1)
        X2 = gf.mul(S1, Y2)
        counters.f2m += 2
        if X1 == 0 or X2 == 0:
            return self._drop(r_full), False
        p1 = int(gf.LOG[X1]) % n
        p2 = int(gf.LOG[X2]) % n
        if p1 == p2:
            return self._drop(r_full), False
        c = r_full.copy()
        c[p1] ^= 1
        c[p2] ^= 1
        return self._drop(c), True
