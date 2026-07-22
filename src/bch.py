"""BCH code object + BPSK/AWGN channel + Berlekamp–Massey hard-decision decoder.

We use a *systematic* BCH encoding by taking a message polynomial f(x) and
producing c(x) = f(x) * g(x). Because the paper's OSD math works with the
generator matrix directly (any full-rank G suffices — the OSD sorts and
column-permutes it anyway), a non-systematic encoding is fine here: what
matters is that our codebook is exactly C_BCH.
"""
from __future__ import annotations

import numpy as np

from .gf import GF, bch_generator_poly, bch_dimension


class BCHCode:
    """Primitive narrow-sense binary BCH code of length n = 2^m - 1."""

    def __init__(self, m: int, t: int):
        self.m = m
        self.t = t
        self.gf = GF(m)
        self.n = self.gf.n
        self.g_poly = bch_generator_poly(self.gf, t)          # coeffs in {0,1}
        self.k = bch_dimension(self.g_poly, self.n)
        self.d_design = 2 * t + 1
        # Build the (k x n) generator matrix. Row i = g(x) * x^i, taken mod 2.
        G = np.zeros((self.k, self.n), dtype=np.int8)
        g = np.array(self.g_poly, dtype=np.int8)
        assert g.size == self.n - self.k + 1
        for i in range(self.k):
            G[i, i:i + g.size] = g
        self.G = G
        # Build a compatible parity-check matrix H (m*t x n) using alpha^{1..2t}.
        # H[i-1, j] = alpha^{i*j} for i=1..2t, j=0..n-1. Each GF(2^m) element is
        # expanded to m binary rows, giving an m*t x n binary H.
        H_ext = np.zeros((self.m * (2 * t), self.n), dtype=np.int64)
        for j in range(self.n):
            for i in range(1, 2 * t + 1):
                e = int(self.gf.pow(2, (i * j) % self.n))  # alpha^(ij)
                # Expand e into m bits (LSB first).
                for b in range(self.m):
                    H_ext[(i - 1) * self.m + b, j] = (e >> b) & 1
        # Reduce H_ext rank down to (n-k) rows via row reduction over GF(2).
        self.H = self._row_reduce_binary(H_ext)
        # Sanity check: G @ H^T = 0 mod 2.
        assert np.all((self.G @ self.H.T) % 2 == 0), "G H^T != 0"

    # ------------------------------------------------------------------
    @staticmethod
    def _row_reduce_binary(M: np.ndarray) -> np.ndarray:
        """Return the (n-k) x n reduced-row-echelon parity-check rows of M."""
        A = M.copy() % 2
        rows, cols = A.shape
        r = 0
        for c in range(cols):
            if r >= rows:
                break
            piv = None
            for i in range(r, rows):
                if A[i, c]:
                    piv = i
                    break
            if piv is None:
                continue
            if piv != r:
                A[[r, piv]] = A[[piv, r]]
            for i in range(rows):
                if i != r and A[i, c]:
                    A[i] ^= A[r]
            r += 1
        # Drop any all-zero rows.
        keep = np.any(A, axis=1)
        return A[keep]

    # ------------------------------------------------------------------
    def encode(self, msg: np.ndarray) -> np.ndarray:
        """Encode a length-k binary vector; returns length-n binary codeword."""
        return (msg @ self.G) % 2

    # ------------------------------------------------------------------
    def bm_decode(self, r_hard: np.ndarray):
        """Berlekamp–Massey hard-decision decoding.

        Returns (decoded_codeword, ok) where ok=False means uncorrectable.
        """
        gf = self.gf
        # Compute 2t syndromes S_i = r(alpha^i).
        S = np.zeros(2 * self.t + 1, dtype=np.int64)
        # r(alpha^i) = sum_j r_j * (alpha^i)^j = sum_j r_j * alpha^{ij}
        idx = np.nonzero(r_hard)[0]
        for i in range(1, 2 * self.t + 1):
            s = 0
            for j in idx:
                s ^= int(gf.pow(2, (i * j) % gf.n))
            S[i] = s
        if not S[1:].any():
            return r_hard.copy(), True

        # Berlekamp–Massey to find error locator Lambda(x).
        L = 0                # current length
        Lam = [1]            # current Lambda
        B = [1]              # previous B(x)
        b = 1                # previous discrepancy
        x_pow_B = [0, 1]     # x * B initial: kept implicit via shift
        m_shift = 1
        for n in range(1, 2 * self.t + 1):
            # Discrepancy: delta = S[n] + sum_{i=1..L} Lam[i] * S[n-i]
            delta = int(S[n])
            for i in range(1, L + 1):
                if i < len(Lam) and Lam[i]:
                    delta ^= gf.mul(int(Lam[i]), int(S[n - i]))
            if delta == 0:
                m_shift += 1
            else:
                # T = Lam - (delta / b) * x^m * B
                coef = gf.div(delta, b)
                # x^m * B
                xmB = [0] * m_shift + list(B)
                # Ensure Lam long enough
                new_len = max(len(Lam), len(xmB))
                T = list(Lam) + [0] * (new_len - len(Lam))
                xmB = xmB + [0] * (new_len - len(xmB))
                for i in range(new_len):
                    T[i] ^= gf.mul(coef, xmB[i])
                if 2 * L <= n - 1:
                    L_new = n - L
                    B = Lam
                    b = delta
                    Lam = T
                    L = L_new
                    m_shift = 1
                else:
                    Lam = T
                    m_shift += 1
        # Chien search: alpha^{-i} is a root of Lam iff position i is an error.
        err_positions = []
        for i in range(self.n):
            # evaluate Lam(alpha^{-i}) = Lam(alpha^{n-i})
            val = 0
            for j, cj in enumerate(Lam):
                if cj == 0:
                    continue
                val ^= gf.mul(int(cj), int(gf.pow(2, ((self.n - i) * j) % gf.n)))
            if val == 0:
                err_positions.append(i)
        # Sanity: number of roots equals deg(Lam)
        if len(err_positions) != L or L > self.t:
            return r_hard.copy(), False
        c = r_hard.copy()
        for p in err_positions:
            c[p] ^= 1
        return c, True


# --------------------------------------------------------------------------
# Channel utilities
# --------------------------------------------------------------------------
def bpsk_modulate(c: np.ndarray) -> np.ndarray:
    """0 -> +1, 1 -> -1."""
    return 1.0 - 2.0 * c.astype(np.float64)


def awgn(x: np.ndarray, sigma: float, rng: np.random.Generator) -> np.ndarray:
    return x + sigma * rng.standard_normal(x.shape)


def sigma_from_ebn0(ebn0_db: float, rate: float) -> float:
    """N0/2 relative to Eb=1 BPSK. sigma^2 = 1/(2 R Eb/N0)."""
    ebn0 = 10 ** (ebn0_db / 10.0)
    return float(np.sqrt(1.0 / (2.0 * rate * ebn0)))


def llr_from_y(y: np.ndarray, sigma: float) -> np.ndarray:
    """LLR log P(y|c=0)/P(y|c=1) for BPSK 0->+1, 1->-1: LLR = 2y/sigma^2."""
    return 2.0 * y / (sigma * sigma)
