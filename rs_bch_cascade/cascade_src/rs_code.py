"""Reed-Solomon code over GF(2^m) with cyclic BCH-like construction.

We use the classical narrow-sense cyclic RS code:
    Generator polynomial: g(x) = prod_{i=1..2t} (x - α^i)
    Codewords: c(x) = m(x) * g(x) (nonsystematic) or systematic form.
    Property: c(α^i) = 0 for i = 1, 2, ..., 2t.

This form makes BM/Chien/Forney straightforward.

For SYSTEMATIC form: c(x) = x^{n-k} * m(x) + [x^{n-k} * m(x) mod g(x)]
- First k symbols = message
- Last n-k symbols = parity (systematic redundancy)

Lagrange interpolation is optional and used for the "Lagrange shared cache"
path (see lagrange_cache.py).

Decoders:
- BM: standard Berlekamp-Massey + Chien search + Forney.
- LCC-BR: Chase decoding on η LRPs.
"""
from __future__ import annotations

import numpy as np
from dataclasses import dataclass

from .upstream import GF, OpCounters


class RSCode:
    """Narrow-sense primitive Reed-Solomon code over GF(2^m).

    Parameters
    ----------
    m : int
        GF exponent (n = 2^m - 1).
    k : int
        Message length in symbols.

    Properties
    ----------
    n : int = 2^m - 1
    t : int = (n - k) // 2, error correction capability.
    d : int = n - k + 1, minimum distance.
    g_poly : list[int] — generator polynomial coefficients (in GF(2^m)).
    """

    def __init__(self, m: int, k: int):
        self.gf = GF(m)
        self.m = m
        self.n = self.gf.n  # 2^m - 1
        self.k = k
        self.d = self.n - self.k + 1
        self.t = (self.d - 1) // 2

        # Locators α^0, α^1, ..., α^{n-1}
        self.alpha_pow = np.array(
            [int(self.gf.EXP[j]) for j in range(self.n)],
            dtype=np.int64)

        # Generator polynomial: g(x) = prod_{i=1..2t} (x - α^i)
        # In char 2 (over GF(2^m)), minus = plus.
        self.g_poly = [1]
        for i in range(1, 2 * self.t + 1):
            root = int(self.gf.EXP[i])
            # Multiply g(x) by (x - root) = (x + root)  in char 2
            self.g_poly = self.gf.poly_mul(self.g_poly, [root, 1])
        assert len(self.g_poly) - 1 == self.n - self.k, \
            f"g_poly degree {len(self.g_poly)-1} != n-k={self.n-self.k}"

    # -------------------------------------------------
    # Encoding
    # -------------------------------------------------
    def encode_systematic(self, msg: np.ndarray) -> np.ndarray:
        """Systematic encoding via polynomial division.

        c(x) = x^{n-k} * m(x) - [x^{n-k} * m(x) mod g(x)]
             = x^{n-k} * m(x) + parity(x)  (char 2)

        Layout: c = [parity_0, ..., parity_{n-k-1}, msg_0, ..., msg_{k-1}]
        (first n-k are parity, last k are message).
        """
        assert msg.size == self.k
        gf = self.gf
        # msg polynomial (deg <= k-1). shift up by (n-k).
        # dividend = msg * x^{n-k}
        # in list form: [0, 0, ..., 0, msg_0, msg_1, ..., msg_{k-1}]
        dividend = [0] * (self.n - self.k) + [int(x) for x in msg]
        # polynomial divison over GF(2^m)
        _quot, rem = gf.poly_divmod(dividend, self.g_poly)
        parity = list(rem) + [0] * (self.n - self.k - len(rem))
        assert len(parity) == self.n - self.k

        c = np.zeros(self.n, dtype=np.int64)
        # c[0..n-k-1] = parity, c[n-k..n-1] = msg
        c[:self.n - self.k] = parity[:self.n - self.k]
        c[self.n - self.k:] = msg
        return c

    def extract_message(self, c: np.ndarray) -> np.ndarray:
        """Extract the message part from a systematic codeword."""
        return c[self.n - self.k:].copy()

    # -------------------------------------------------
    # Berlekamp-Massey hard-decision decoder
    # -------------------------------------------------
    def bm_decode(self, r: np.ndarray, counters: OpCounters = None):
        """Standard RS BM + Chien + Forney.

        r: received word (length n), symbols in GF(2^m).
        Returns (decoded_codeword, ok).
        """
        gf = self.gf
        n = self.n
        t = self.t
        if counters is None:
            counters = OpCounters()

        # 1. Compute 2t syndromes S_i = r(α^i) for i = 1..2t
        # Vectorized: for each i, we need sum_j r[j] * α^{i*j} over nonzero r[j]
        S = np.zeros(2 * t + 1, dtype=np.int64)
        nz = np.nonzero(r)[0]
        if nz.size == 0:
            return r.copy(), True
        # For each i, compute the syndrome as XOR over nonzero r contributions
        for i in range(1, 2 * t + 1):
            # α^{i*j} for j in nz
            exps = (i * nz) % gf.n
            alpha_ij = gf.EXP[exps]
            # r[j] * α^{i*j} — element-wise GF mult
            # Use LOG table for speed
            r_nz = r[nz].astype(np.int64)
            r_log = gf.LOG[r_nz]
            alpha_log = gf.LOG[alpha_ij]
            prod_log = (r_log + alpha_log) % gf.n
            prod = gf.EXP[prod_log]
            S[i] = int(np.bitwise_xor.reduce(prod))
        counters.f2m += int(nz.size) * 2 * t

        if not S[1:].any():
            return r.copy(), True

        # 2. BM iteration
        L = 0
        Lam = [1]
        B = [1]
        b = 1
        m_shift = 1
        for k_iter in range(1, 2 * t + 1):
            delta = int(S[k_iter])
            for i in range(1, L + 1):
                if i < len(Lam) and Lam[i]:
                    delta ^= gf.mul(int(Lam[i]), int(S[k_iter - i]))
                    counters.f2m += 1
            if delta == 0:
                m_shift += 1
            else:
                coef = gf.div(delta, b)
                counters.f2m += 1
                xmB = [0] * m_shift + list(B)
                new_len = max(len(Lam), len(xmB))
                T = list(Lam) + [0] * (new_len - len(Lam))
                xmB = xmB + [0] * (new_len - len(xmB))
                for i in range(new_len):
                    T[i] ^= gf.mul(coef, xmB[i])
                    counters.f2m += 1
                if 2 * L <= k_iter - 1:
                    L_new = k_iter - L
                    B = Lam
                    b = delta
                    Lam = T
                    L = L_new
                    m_shift = 1
                else:
                    Lam = T
                    m_shift += 1

        # 3. Chien search: find roots of Λ(x) = α^{-p} for each error position p
        # Vectorized: for each position i in 0..n-1, evaluate Λ(α^{n-i})
        Lam_arr = np.array(Lam, dtype=np.int64)
        Lam_nz_deg = np.where(Lam_arr != 0)[0]
        if Lam_nz_deg.size == 0:
            return r.copy(), False
        err_positions = []
        for i in range(n):
            # Evaluate Λ(α^{n-i}) = sum_{d} Lam[d] * α^{d*(n-i)}
            exps = ((n - i) * Lam_nz_deg) % gf.n
            log_lam = gf.LOG[Lam_arr[Lam_nz_deg]]
            log_alpha = gf.LOG[gf.EXP[exps]]
            prod_log = (log_lam + log_alpha) % gf.n
            prod = gf.EXP[prod_log]
            val = int(np.bitwise_xor.reduce(prod))
            if val == 0:
                err_positions.append(i)
        counters.f2m += n * int(Lam_nz_deg.size)

        if len(err_positions) != L or L > t:
            return r.copy(), False

        # 4. Forney error evaluation
        # Ω(x) = [S(x) * Λ(x)] mod x^{2t+1}
        Sx = [0] + [int(S[i]) for i in range(1, 2 * t + 1)]
        Omega = [0] * (2 * t + 1)
        for i in range(len(Sx)):
            if Sx[i] == 0:
                continue
            for j in range(len(Lam)):
                if i + j <= 2 * t and Lam[j]:
                    Omega[i + j] ^= gf.mul(Sx[i], Lam[j])
                    counters.f2m += 1

        # Λ'(x) — formal derivative in char 2 (odd-degree terms only)
        Lam_prime = [0] * len(Lam)
        for i in range(1, len(Lam)):
            if i % 2 == 1:
                Lam_prime[i - 1] = Lam[i]

        c = r.copy()
        for p in err_positions:
            x_inv = int(gf.EXP[(n - p) % gf.n])
            omega_val = 0
            for i, oi in enumerate(Omega):
                if oi == 0:
                    continue
                omega_val ^= gf.mul(oi, int(gf.EXP[(i * (n - p)) % gf.n]))
                counters.f2m += 1
            lam_prime_val = 0
            for i, li in enumerate(Lam_prime):
                if li == 0:
                    continue
                lam_prime_val ^= gf.mul(li, int(gf.EXP[(i * (n - p)) % gf.n]))
                counters.f2m += 1
            if lam_prime_val == 0:
                return r.copy(), False
            # For narrow-sense: e_p = -α^p * Ω(α^{-p}) / Λ'(α^{-p})
            # In char 2, - = +. So e_p = α^p * Ω(x_inv) / Λ'(x_inv)
            err_val = gf.div(omega_val, lam_prime_val)
            # For narrow-sense RS: error value directly (without α^p factor
            # when using this specific formulation). Let's verify.
            # Actually the Forney formula with S_i = r(α^i) starting at i=1:
            #   e_p = -Ω(α^{-p}) / [α^{-p} * Λ'(α^{-p})]
            # In char 2 that simplifies to Ω(x_inv) / (x_inv * Λ'(x_inv))
            # which equals Ω(x_inv) * α^p / Λ'(x_inv).
            alpha_p = int(gf.EXP[p % gf.n])
            err_val = gf.mul(alpha_p, err_val)
            counters.f2m += 1
            c[p] ^= err_val

        return c, True

    # -------------------------------------------------
    # LCC-BR soft decoder (simplified via BM on Chase test-vectors)
    # -------------------------------------------------
    def lcc_br_decode(self, r_soft: np.ndarray, reliability: np.ndarray,
                       eta: int, counters: OpCounters = None):
        """Chase-BM style soft decoder.

        r_soft: n received symbols (over GF(2^m)).
        reliability: per-symbol reliability score (higher = more reliable).
        eta: number of least-reliable positions used for 2^η test-vectors.

        Simplified LCC-BR: for each test-vector, run BM on modified r.
        Pick the best (weighted matching score with reliability).
        """
        from itertools import product

        n = self.n
        t = self.t
        gf = self.gf
        if counters is None:
            counters = OpCounters()

        lrp_indices = np.argsort(reliability)[:eta]

        best_c = None
        best_score = -np.inf
        n_tvs = 0

        for bits in product([0, 1], repeat=eta):
            n_tvs += 1
            r_test = r_soft.copy()
            for i, b in enumerate(bits):
                if b:
                    r_test[lrp_indices[i]] ^= 1
            c_hat, ok = self.bm_decode(r_test, counters)
            if not ok:
                continue
            match_mask = (c_hat == r_soft)
            score = float(np.sum(reliability[match_mask]))
            if score > best_score:
                best_score = score
                best_c = c_hat

        if best_c is None:
            c_hat, ok = self.bm_decode(r_soft, counters)
            return c_hat if ok else r_soft, ok

        if not hasattr(counters, "n_tvs"):
            counters.n_tvs = 0
        counters.n_tvs += n_tvs

        return best_c, True
