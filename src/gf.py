"""GF(2^m) finite-field arithmetic used throughout the LLOSD reproduction.

Elements are represented as integers 0..2^m-1 with the natural binary
polynomial representation (bit i is coefficient of alpha^i where alpha is a
root of the primitive polynomial).

We build LOG/EXP tables at construction; multiplication/inverse/pow are then
O(1). Addition is XOR.
"""
from __future__ import annotations

import numpy as np


# Primitive polynomials for common m, as integers (bit i = coeff of x^i).
# The extra leading bit at position m encodes the x^m term.
PRIM_POLY = {
    1: 0b11,             # x + 1
    2: 0b111,            # x^2 + x + 1
    3: 0b1011,           # x^3 + x + 1
    4: 0b10011,          # x^4 + x + 1
    5: 0b100101,         # x^5 + x^2 + 1
    6: 0b1000011,        # x^6 + x + 1
    7: 0b10001001,       # x^7 + x^3 + 1
    8: 0b100011101,      # x^8 + x^4 + x^3 + x^2 + 1
}


class GF:
    """Finite field GF(2^m) built from an irreducible primitive polynomial."""

    def __init__(self, m: int):
        if m not in PRIM_POLY:
            raise ValueError(f"No primitive polynomial for m={m}")
        self.m = m
        self.n = (1 << m) - 1  # field order minus one
        self.prim = PRIM_POLY[m]
        # Build EXP/LOG tables.
        self.EXP = np.zeros(2 * self.n + 2, dtype=np.int64)
        self.LOG = np.full(1 << m, -1, dtype=np.int64)
        x = 1
        for i in range(self.n):
            self.EXP[i] = x
            self.LOG[x] = i
            x <<= 1
            if x & (1 << m):
                x ^= self.prim
        # Duplicate for mod-free lookup: EXP[i + n] == EXP[i].
        for i in range(self.n, 2 * self.n + 2):
            self.EXP[i] = self.EXP[i - self.n]
        # By convention: LOG[0] is -inf; we use -1 as sentinel and guard callers.

    # --- scalar ops --------------------------------------------------------
    def add(self, a: int, b: int) -> int:
        return a ^ b

    def mul(self, a: int, b: int) -> int:
        if a == 0 or b == 0:
            return 0
        return int(self.EXP[self.LOG[a] + self.LOG[b]])

    def inv(self, a: int) -> int:
        if a == 0:
            raise ZeroDivisionError("inverse of zero in GF(2^m)")
        return int(self.EXP[self.n - self.LOG[a]])

    def div(self, a: int, b: int) -> int:
        if a == 0:
            return 0
        if b == 0:
            raise ZeroDivisionError("divide by zero in GF(2^m)")
        return int(self.EXP[self.LOG[a] - self.LOG[b] + self.n])

    def pow(self, a: int, e: int) -> int:
        if a == 0:
            return 0 if e > 0 else 1
        return int(self.EXP[(self.LOG[a] * e) % self.n])

    # --- vector ops (numpy) ------------------------------------------------
    def vadd(self, A: np.ndarray, B: np.ndarray) -> np.ndarray:
        return np.bitwise_xor(A, B)

    def vmul(self, A: np.ndarray, B: np.ndarray) -> np.ndarray:
        """Element-wise GF multiplication using EXP/LOG tables."""
        A = np.asarray(A, dtype=np.int64)
        B = np.asarray(B, dtype=np.int64)
        out = np.zeros(np.broadcast(A, B).shape, dtype=np.int64)
        mask = (A != 0) & (B != 0)
        if mask.any():
            la = self.LOG[A[mask]]
            lb = self.LOG[B[mask]]
            out[mask] = self.EXP[la + lb]
        return out

    def vinv(self, A: np.ndarray) -> np.ndarray:
        A = np.asarray(A, dtype=np.int64)
        if np.any(A == 0):
            raise ZeroDivisionError("inverse of zero in vinv")
        return self.EXP[self.n - self.LOG[A]].astype(np.int64)

    # --- polynomial ops over GF(2^m) --------------------------------------
    def poly_eval(self, coeffs, x: int) -> int:
        """Horner evaluate p(x) where coeffs[i] is the coefficient of x^i."""
        y = 0
        for c in reversed(list(coeffs)):
            y = self.mul(y, x) ^ int(c)
        return y

    def poly_mul(self, a, b):
        a = list(a); b = list(b)
        out = [0] * (len(a) + len(b) - 1)
        for i, ai in enumerate(a):
            if ai == 0:
                continue
            for j, bj in enumerate(b):
                if bj == 0:
                    continue
                out[i + j] ^= self.mul(ai, bj)
        return out

    def poly_add(self, a, b):
        la, lb = len(a), len(b)
        L = max(la, lb)
        out = [0] * L
        for i in range(L):
            va = a[i] if i < la else 0
            vb = b[i] if i < lb else 0
            out[i] = va ^ vb
        return out

    def poly_divmod(self, num, den):
        """Return (quotient, remainder) of poly division over GF(2^m)."""
        num = list(num)
        den = list(den)
        while den and den[-1] == 0:
            den.pop()
        if not den:
            raise ZeroDivisionError("zero divisor")
        lead_inv = self.inv(den[-1])
        quot = [0] * max(0, len(num) - len(den) + 1)
        rem = list(num)
        for i in range(len(quot) - 1, -1, -1):
            if len(rem) <= i + len(den) - 1:
                continue
            coeff = self.mul(rem[i + len(den) - 1], lead_inv)
            quot[i] = coeff
            for j in range(len(den)):
                rem[i + j] ^= self.mul(coeff, den[j])
        while len(rem) > 0 and rem[-1] == 0:
            rem.pop()
        return quot, rem


def bch_generator_poly(gf: GF, t: int):
    """Return the narrow-sense binary primitive BCH generator polynomial g(x).

    g(x) is the least-degree binary polynomial that has alpha, alpha^2, ...,
    alpha^{2t} as roots. Because g is over GF(2), it is the product of the
    minimal polynomials of the odd powers alpha^1, alpha^3, ..., alpha^{2t-1}
    (each covers its whole 2-conjugate cyclotomic class).
    """
    seen = set()
    g = [1]  # constant 1
    for i in range(1, 2 * t + 1):
        # cyclotomic coset of i modulo n
        coset = set()
        j = i
        while j not in coset:
            coset.add(j)
            j = (j * 2) % gf.n
        rep = min(coset)
        if rep in seen:
            continue
        seen.add(rep)
        # Minimal polynomial: prod_{s in coset} (x - alpha^s).
        # Since we are in char 2, minus becomes plus.
        m_poly = [1]
        for s in coset:
            root = int(gf.EXP[s])
            m_poly = gf.poly_mul(m_poly, [root, 1])
        # m_poly should be over GF(2); its coefficients happen to lie in {0,1}.
        # Convert to plain 0/1 ints.
        m_poly_bin = [int(c) & 1 for c in m_poly]
        # But actually all coefficients ARE 0 or 1; assert as a sanity check.
        for c in m_poly:
            assert c in (0, 1), f"minimal poly not binary: {m_poly}"
        g = gf.poly_mul(g, m_poly_bin)
        # Reduce back to binary (poly_mul is XOR-based already, but safe cast):
        g = [int(c) & 1 for c in g]
    return g  # binary polynomial, deg = n - k


def bch_dimension(g_poly, n: int) -> int:
    """k = n - deg(g). Assumes g_poly is binary (coefficients 0/1)."""
    deg = len(g_poly) - 1
    while deg > 0 and g_poly[deg] == 0:
        deg -= 1
    return n - deg
