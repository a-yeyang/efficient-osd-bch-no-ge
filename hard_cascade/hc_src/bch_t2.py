"""BCH(n, k) t=2 hard-decision decoders.

Two decoders implemented per Lagendijk et al. 2026:
  1. Conventional: Berlekamp-Massey (simplified for t≤4 binary BCH) + Chien search
  2. Direct root finding: closed-form solution using LUT (Table I of paper)

Both operate on hard-decision received words and produce corrected codewords.

For binary primitive narrow-sense BCH(n, k, t) with n = 2^m - 1:
  Only ODD syndromes S_1, S_3, ..., S_{2t-1} are needed (because S_{2i} = S_i^2 in GF(2^m)).

For t = 2 the ELP is:
    Λ(x) = 1 + Λ_1 x + Λ_2 x^2
where by BM:
    Λ_1 = S_1
    Λ_2 = (S_1^3 + S_3) / S_1     if S_1 != 0

Error scenarios (Lagendijk Table I):
    0 errors:  S_1 = 0 AND S_3 = 0
    1 error:   S_1 != 0 AND S_1^3 = S_3         (only need Λ_1)
    2 errors:  S_1 != 0 AND S_1^3 != S_3        (full quadratic)
    3+ errors: S_1 = 0 AND S_3 != 0             (uncorrectable)
"""
from __future__ import annotations

import numpy as np
from dataclasses import dataclass

from .upstream import GF, OpCounters


class BCHt2Code:
    """Binary primitive narrow-sense BCH(n, k, t=2) code with two hard decoders."""

    def __init__(self, m: int):
        self.m = m
        self.gf = GF(m)
        self.n = self.gf.n              # 2^m - 1
        self.t = 2
        self.d_design = 2 * self.t + 1  # 5

        # Generator polynomial: g(x) = m_1(x) * m_3(x)
        # where m_i is the minimal polynomial of α^i over GF(2).
        # For m=8: m_1 has degree 8 (cyclotomic class of 1 mod 255 has 8 elems);
        # m_3 has degree 8. So g_deg = 16, k = n-16.
        seen = set()
        g_poly = [1]
        for i in [1, 3]:
            coset = self._cyclotomic_coset(i)
            rep = min(coset)
            if rep in seen: continue
            seen.add(rep)
            m_poly = [1]
            for s in coset:
                root = int(self.gf.EXP[s])
                m_poly = self.gf.poly_mul(m_poly, [root, 1])
            m_poly_bin = [int(c) & 1 for c in m_poly]
            g_poly = self.gf.poly_mul(g_poly, m_poly_bin)
            g_poly = [int(c) & 1 for c in g_poly]
        self.g_poly = g_poly
        self.k = self.n - (len(g_poly) - 1)
        assert self.k == self.n - 2 * self.m, \
            f"k={self.k} should be n-2m={self.n - 2*self.m} for BCH t=2"

        # Systematic G matrix (k x n): row i encodes msg unit vector e_i
        # Layout: c(x) = x^{n-k} * m(x) + [x^{n-k} m(x) mod g(x)]
        # c = [parity(n-k bits), msg(k bits)]
        self.G = self._build_G()
        # Parity check matrix H (from BCH structure)
        # Not needed for the two decoders below

        # Precompute LUT for Direct t=2 decoder
        self.lut_A = self._build_lut_A()

    def _cyclotomic_coset(self, i):
        coset = set()
        j = i
        while j not in coset:
            coset.add(j)
            j = (j * 2) % self.gf.n
        return coset

    def _build_G(self):
        """Build systematic G: k x n binary matrix."""
        n, k = self.n, self.k
        n_minus_k = n - k
        G = np.zeros((k, n), dtype=np.int8)
        # For each message unit vector e_i, encode as c = [parity_i, e_i]
        # where parity_i = [e_i * x^{n-k} mod g(x)]
        gf = self.gf
        g = list(self.g_poly)
        for i in range(k):
            # m(x) = x^i, shift by n-k -> polynomial with '1' at position n-k+i
            dividend = [0] * (n_minus_k + i) + [1] + [0] * (k - i - 1)
            # Trim to length n
            dividend = dividend[:n]
            _q, rem = gf.poly_divmod(dividend, g)
            parity = list(rem) + [0] * (n_minus_k - len(rem))
            G[i, :n_minus_k] = parity[:n_minus_k]
            G[i, n_minus_k + i] = 1
        return G

    def encode(self, msg: np.ndarray) -> np.ndarray:
        """Encode a k-bit message → n-bit codeword. c = msg @ G mod 2."""
        assert msg.size == self.k
        return (msg.astype(np.int64) @ self.G) % 2

    def extract_message(self, cw: np.ndarray) -> np.ndarray:
        """Extract message bits from systematic codeword."""
        n_minus_k = self.n - self.k
        return cw[n_minus_k:]

    # =======================================================================
    # Syndromes (shared by both decoders)
    # =======================================================================
    def compute_syndromes(self, r: np.ndarray, counters: OpCounters = None):
        """Compute S_1 and S_3 (odd syndromes; S_2, S_4 not needed for binary BCH).

        Returns (S1, S3) as GF(2^m) elements (int).
        """
        if counters is None:
            counters = OpCounters()
        gf = self.gf
        n = self.n
        nz = np.nonzero(r)[0]
        if nz.size == 0:
            return 0, 0
        # S_1 = sum_{j: r[j]=1} α^j
        s1 = int(np.bitwise_xor.reduce(gf.EXP[nz % gf.n]))
        # S_3 = sum_{j: r[j]=1} α^{3j}
        s3 = int(np.bitwise_xor.reduce(gf.EXP[(3 * nz) % gf.n]))
        # Ops accounting: 2n XOR of m-bit values → 2n F_{2^m} adds
        counters.f2m += 2 * n
        return s1, s3

    # =======================================================================
    # Decoder 1: Conventional (BM + Chien)
    # =======================================================================
    def decode_conventional(self, r: np.ndarray, counters: OpCounters = None):
        """Berlekamp-Massey (simplified for t=2 binary BCH) + Chien search.

        Returns (decoded_codeword, ok).
        """
        gf = self.gf
        n = self.n
        if counters is None:
            counters = OpCounters()

        S1, S3 = self.compute_syndromes(r, counters)

        # BM for t=2 binary BCH: 4 syndromes (S1, S2=S1², S3, S4=S3² ... no wait).
        # Actually for binary BCH we only need S_1 and S_3 because
        # S_{2i} = S_i^2, and the BM operates on the 2t=4 syndrome sequence
        # (S_1, S_2, S_3, S_4) but only S_1, S_3 are informative.
        #
        # From Lagendijk paper (Table I):
        #   0 errors: S_1 = 0 and S_3 = 0
        #   1 error:  S_1 != 0 and S_1^3 = S_3
        #   2 errors: S_1 != 0 and S_1^3 != S_3
        #   3+:       S_1 = 0 and S_3 != 0
        S1_cubed = gf.mul(gf.mul(S1, S1), S1)  # 2 mults
        counters.f2m += 2

        if S1 == 0 and S3 == 0:
            return r.copy(), True

        if S1 == 0 and S3 != 0:
            # Uncorrectable (odd errors, but S_1 = 0 means even parity... failure)
            return r.copy(), False

        # S_1 != 0
        if S1_cubed == S3:
            # 1 error case: For 1 error at position p, S_1 = α^p, S_3 = α^{3p} = S_1^3.
            # Error position: p = log_α(S_1)
            log_s1 = int(gf.LOG[S1])
            p = log_s1 % n
            c = r.copy()
            c[p] ^= 1
            return c, True

        # 2 error case: monic ELP
        # Λ_monic(x) = x^2 + S_1 x + (S_1^3 + S_3)/S_1
        # (This is the form used by Lagendijk direct root finding.)
        # Roots of Λ_monic are α^{p_1}, α^{p_2} (elementary symmetric of
        # error locators X_i = α^{p_i}, where p_i are the error positions).
        # Verify: Λ_monic(x) = (x - X_1)(x - X_2) = x^2 - (X_1+X_2)x + X_1 X_2
        # in char 2: x^2 + (X_1+X_2)x + X_1 X_2
        # S_1 = X_1 + X_2 (elementary symmetric power sum p_1 = e_1)
        # Constant term = X_1 X_2 = e_2 = (S_1^3 + S_3)/S_1 (Newton's identity)
        c0 = gf.div(gf.add(S1_cubed, S3), S1)  # constant term of Λ_monic
        counters.f2m += 2

        # Chien-style search on monic form:
        # For each candidate position i in 0..n-1, check if X = α^i is a root:
        #   Λ_monic(α^i) = α^{2i} + S_1 α^i + c0
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
            return r.copy(), False

        c = r.copy()
        for p in err_positions:
            c[p] ^= 1
        return c, True

    # =======================================================================
    # Decoder 2: Direct root finding (Lagendijk §III-A)
    # =======================================================================
    def _build_lut_A(self):
        """Build the {}_A LUT (Lagendijk Table I).

        The LUT maps k → the two roots of A(X) = X^2 + X + k.
        We build it by iterating over all X and computing k = X^2 + X for each,
        then storing the (X, X+1) pair.

        LUT format:
          - Entry index by k (0 to 2^m - 1)
          - Each entry: (X_root1, X_root2) — but only half of k values have roots.
          - Roots exist iff Tr(k) = 0 (trace over GF(2^m) → GF(2))
          - When roots exist, the two roots are X and X+1.
        """
        gf = self.gf
        n_field = 1 << self.m  # 2^m
        lut = np.zeros((n_field, 2), dtype=np.int64)
        lut_valid = np.zeros(n_field, dtype=bool)

        for X in range(n_field):
            X2 = gf.mul(X, X)
            k = X2 ^ X  # X^2 + X = X^2 XOR X in char 2
            if not lut_valid[k]:
                lut[k, 0] = X
                lut[k, 1] = X ^ 1  # X + 1 in char 2
                lut_valid[k] = True

        self.lut_A_valid = lut_valid
        return lut

    def decode_direct(self, r: np.ndarray, counters: OpCounters = None):
        """Direct root finding for t=2 BCH (Lagendijk §III-A).

        Complete workflow:
          1. Compute S_1, S_3 (same as conventional)
          2. Determine error count from (S_1, S_3) conditions
          3. For 2-error case: transform Λ into A(Y) = Y^2 + Y + k form,
             lookup roots in LUT, then invert.

        Returns (decoded_codeword, ok).
        """
        gf = self.gf
        n = self.n
        if counters is None:
            counters = OpCounters()

        S1, S3 = self.compute_syndromes(r, counters)
        S1_cubed = gf.mul(gf.mul(S1, S1), S1)
        counters.f2m += 2

        if S1 == 0 and S3 == 0:
            return r.copy(), True
        if S1 == 0 and S3 != 0:
            return r.copy(), False

        if S1_cubed == S3:
            # 1 error: For 1 error at position p, S_1 = α^p → p = log_α(S_1)
            log_s1 = int(gf.LOG[S1])
            p = log_s1 % n
            c = r.copy()
            c[p] ^= 1
            return c, True

        # 2 errors: monic ELP
        # Λ_monic(x) = x^2 + S_1 x + (S_1^3 + S_3)/S_1
        # Substitute x = S_1 * Y: gives S_1^2 (Y^2 + Y) + (S_1^3 + S_3)/S_1
        # Divide by S_1^2:  A(Y) = Y^2 + Y + (S_1^3+S_3)/S_1^3
        # LUT gives roots Y_1, Y_2 of A(Y); roots of Λ are X_i = S_1 * Y_i.
        # Since Λ_monic's roots are α^{p_i} (see comment in decode_conventional),
        # error positions p_i = log_α(X_i).
        numerator = gf.add(S1_cubed, S3)  # XOR
        k_lut = gf.div(numerator, S1_cubed)
        counters.f2m += 1

        # LUT lookup: single cycle in hardware
        if not self.lut_A_valid[k_lut]:
            return r.copy(), False

        Y1 = int(self.lut_A[k_lut, 0])
        Y2 = int(self.lut_A[k_lut, 1])

        # Roots of Λ_monic: X_i = S_1 * Y_i
        X1 = gf.mul(S1, Y1)
        X2 = gf.mul(S1, Y2)
        counters.f2m += 2

        if X1 == 0 or X2 == 0:
            return r.copy(), False

        # Error positions: X_i = α^{p_i} → p_i = log_α(X_i)
        p1 = int(gf.LOG[X1]) % n
        p2 = int(gf.LOG[X2]) % n

        if p1 == p2:
            return r.copy(), False

        c = r.copy()
        c[p1] ^= 1
        c[p2] ^= 1
        return c, True

    # =======================================================================
    # Hardware clock cycle model (Lagendijk Table VI extrapolated)
    # =======================================================================
    def cycles_conventional(self) -> int:
        """Latency in clock cycles for conventional BM + Chien search.

        Per Lagendijk Table VI (n=256, t=2): 8 cycles.
        Breakdown:
          - Syndrome compute: 1 cycle (n-parallel XOR reduce)
          - BM iteration: 2t = 4 cycles (fully parallel BM)
          - Chien search: 1 cycle (n-way parallel evaluation)
          - Error correction + margin: 2 cycles
        We use paper's number directly for n=256; for other n we assume
        the same structural cycle count.
        """
        return 2 * self.t + 2 + 2  # 8 for t=2

    def cycles_direct(self) -> int:
        """Latency for Direct root finding.

        Per Lagendijk Table VI (n=256, t=2): 3 cycles.
        Breakdown:
          - Syndrome compute: 1 cycle
          - Precomputation (k_lut): 1 cycle
          - LUT lookup + inversion + correction: 1 cycle
        """
        return 3


class BCHt2CodeN128(BCHt2Code):
    """Alias for BCH(127, k, t=2) — n=127 in GF(2^7)."""
    def __init__(self):
        super().__init__(m=7)
