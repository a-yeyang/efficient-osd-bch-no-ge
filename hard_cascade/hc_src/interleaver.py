"""Matrix interleaving for RS+BCH concatenated codes (advisor's 4th task).

Implements the interleaving scheme of Liu/Song/Wang, "A Novel Interleaving
Scheme for Concatenated Codes on Burst-Error Channel," APCC 2022, applied to our
Config 1 (KP4): the paper's *exact* code — x=2 RS(544,514)/GF(2^10) outer
codewords feeding y=80 BCH(144,136,t=1)/GF(2^8) inner codewords.

Two independent permutations, matching the paper's decomposition:

  * Step 1 (RS preprocessing): reorder the x RS codewords' symbols so that the
    residual symbol errors left after BCH decoding are spread *evenly across the
    x RS codewords* (neither codeword alone exceeds its t=15 budget).  We
    interleave the two codewords symbol-by-symbol: A0,B0,A1,B1,...  — a clean,
    invertible realisation of the paper's odd/even part reordering.

  * Steps 2-3 (matrix interleaving + channel order): a block interleaver over the
    y BCH codewords with configurable *column width* w (unit_bits).  The y coded
    codewords are the rows; we transmit unit-by-unit down the rows (column-major
    over units), so a channel burst of B consecutive bits is dispersed across
    ~B/w different BCH codewords — each seeing at most one unit of errors.

Column-width settings from the paper:
  * w = 1 bit    : maximal dispersion → best BCH protection on burst channels.
  * w = 2 bits   : one PAM4 signal; adjacent PAM4 symbols land in different BCH
                   codewords, but a symbol's 2 Gray bits stay together (a burst
                   slip usually flips only 1 of them → still ≤ t=1).
  * w = m bits   : one RS symbol (m = m_rs = 10); preserves RS-symbol locality
                   (good for RS) at the cost of less BCH dispersion.

Every scheme here is a *pure bit-position permutation*, so interleave followed by
deinterleave is the identity by construction (unit-tested).  On a memoryless
(AWGN) channel any permutation is error-distribution-neutral, so all schemes
coincide — the gain appears only on burst channels.
"""
from __future__ import annotations

import math
import numpy as np


def build_tx_perm(n_blocks: int, block_bits: int, unit_bits: int) -> np.ndarray:
    """Block-interleaver permutation over `n_blocks` BCH codewords.

    Source layout: codeword b, bit j  ->  source index  b*block_bits + j.
    Transmit order: for each unit-column, walk all codewords (column-major),
    emitting `unit_bits` bits at a time (ragged last column when unit_bits does
    not divide block_bits).

    Returns `perm` with  tx_stream = source[perm]  (perm[tx_pos] = src_idx).
    `unit_bits >= block_bits` yields the identity (no interleaving).
    """
    n_units = math.ceil(block_bits / unit_bits)
    perm = np.empty(n_blocks * block_bits, dtype=np.int64)
    w = 0
    for uc in range(n_units):
        lo = uc * unit_bits
        hi = min(lo + unit_bits, block_bits)
        for b in range(n_blocks):
            base = b * block_bits
            for j in range(lo, hi):
                perm[w] = base + j
                w += 1
    assert w == perm.size
    return perm


class MatrixInterleaver:
    """Transmission-order block interleaver (paper Steps 2-3).

    mode:
      'none'      : sequential transmission (unit = whole codeword) — no spread.
      'existing'  : symbol-level interleave (unit = m_rs bits), no Step-1.
      'w1bit'     : proposed, unit = 1 bit.
      'w2bits'    : proposed, unit = 2 bits (one PAM4 signal).
      'w1sym'     : proposed, unit = m_rs bits (one RS symbol).
    ('existing' vs 'w1sym' share the unit width but differ in Step-1, which is
     applied by the codec, not here.)
    """

    _UNIT = {"none": None, "existing": "m", "w1bit": 1, "w2bits": 2, "w1sym": "m"}

    def __init__(self, n_blocks: int, block_bits: int, m_rs: int, mode: str):
        assert mode in self._UNIT, f"unknown interleaver mode {mode}"
        self.n_blocks = n_blocks
        self.block_bits = block_bits
        self.m_rs = m_rs
        self.mode = mode

        if mode == "none":
            unit = block_bits            # identity
        else:
            u = self._UNIT[mode]
            unit = m_rs if u == "m" else u
        self.unit_bits = unit
        self.perm = build_tx_perm(n_blocks, block_bits, unit)
        # inverse: source[perm] = tx  =>  deinterleave places rx back by perm
        self.n = self.perm.size

    def interleave(self, coded_bits: np.ndarray) -> np.ndarray:
        assert coded_bits.size == self.n
        return coded_bits[self.perm]

    def deinterleave(self, rx_bits: np.ndarray) -> np.ndarray:
        assert rx_bits.size == self.n
        out = np.empty_like(rx_bits)
        out[self.perm] = rx_bits
        return out

    def uses_step1(self) -> bool:
        """Proposed schemes apply the Step-1 RS preprocessing; 'existing'/'none' do not."""
        return self.mode in ("w1bit", "w2bits", "w1sym")


def rs_symbol_interleave(symbol_streams: list[np.ndarray]) -> np.ndarray:
    """Step-1: interleave x RS codewords' symbols round-robin (A0,B0,A1,B1,...).

    All streams must have equal length.  Spreads residual (post-BCH) symbol
    errors evenly across the x RS codewords.
    """
    x = len(symbol_streams)
    L = symbol_streams[0].size
    for s in symbol_streams:
        assert s.size == L
    out = np.empty(x * L, dtype=np.int64)
    for i in range(x):
        out[i::x] = symbol_streams[i]
    return out


def rs_symbol_deinterleave(flat: np.ndarray, x: int) -> list[np.ndarray]:
    """Inverse of rs_symbol_interleave: split back into x symbol streams."""
    return [flat[i::x].copy() for i in range(x)]
