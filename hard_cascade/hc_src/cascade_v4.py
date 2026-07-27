"""v4 cascade codec with matrix interleaving (advisor's 4th task, 2026-07-27).

Adds the Liu/Song/Wang APCC-2022 interleaving scheme to the v3 RS+BCH cascade.
Built on Config 1 (KP4) because that is the paper's *exact* code:

    x = 2  RS(544,514,t=15)/GF(2^10)  outer codewords
    y = 80 BCH(144,136,t=1)/GF(2^8)   inner codewords
    (2*544*10 = 10880 bits = 80*136, an exact fit — no padding.)

Encode chain (interleaved):
    2 msgs -> RS encode each -> [Step-1 symbol interleave] -> serialise to bits
           -> 80 BCH(136->144) encodes -> [Steps 2-3 bit interleave] -> channel
Decode chain: inverse in reverse order.

Modes (interleaver): none / existing / w1bit / w2bits / w1sym.
'none' reproduces a sequential (no-interleave) cascade; on AWGN every mode is
statistically identical (permutation-neutral), so differences appear only on the
burst channel.
"""
from __future__ import annotations

import numpy as np

from .upstream import OpCounters
from .shortened_codes import ShortenedRSCode, BinaryBCH
from .interleaver import (
    MatrixInterleaver, rs_symbol_interleave, rs_symbol_deinterleave,
)
from .burst_channel import run_pam4_channel_burst
from .cascade_v3 import run_pam4_channel_hard


# The paper's frame: 2 RS codewords of Config 1 feeding 80 BCH(144,136).
PAPER_M_RS = 10
PAPER_N_RS = 544
PAPER_K_RS = 514
PAPER_M_BCH = 8
PAPER_T_BCH = 1
PAPER_N_BCH = 144
PAPER_K_BCH = 136
PAPER_X_RS = 2          # outer RS codewords per frame
PAPER_Y_BCH = 80        # inner BCH codewords per frame


# ---------------------------------------------------------------------------
# Fast vectorised RS parity encoder (matches ShortenedRSCode.encode_systematic
# bit-for-bit).  Parity is GF-linear in the message; we precompute the parity
# generator matrix P (k x n_parity) once by encoding unit symbols, then encode
# via a log-table GF matrix product.  Shared across all interleaver modes since
# they use the identical RS code.
# ---------------------------------------------------------------------------
class _FastRSParity:
    def __init__(self, rs: ShortenedRSCode):
        self.rs = rs
        gf = rs.gf
        self.gf = gf
        k, npar = rs.k, rs.n_parity
        P = np.zeros((k, npar), dtype=np.int64)
        e = np.zeros(k, dtype=np.int64)
        for i in range(k):
            e[:] = 0
            e[i] = 1
            P[i] = rs.encode_systematic(e)[:npar]     # unit-symbol parity
        self.P = P
        self.LOGP = np.where(P > 0, gf.LOG[P], -1)    # -1 marks a zero entry
        self.n_field = gf.n

    def encode_systematic(self, msg: np.ndarray) -> np.ndarray:
        gf = self.gf
        msg = msg.astype(np.int64)
        nz = msg != 0
        parity = np.zeros(self.rs.n_parity, dtype=np.int64)
        if nz.any():
            logm = gf.LOG[msg[nz]]                     # (nnz,)
            logp = self.LOGP[nz]                       # (nnz, n_parity)
            valid = logp >= 0
            prod = np.zeros_like(logp)
            s = (logm[:, None] + logp) % self.n_field
            prod[valid] = gf.EXP[s[valid]]
            parity = np.bitwise_xor.reduce(prod, axis=0).astype(np.int64)
        c = np.empty(self.rs.n, dtype=np.int64)
        c[:self.rs.n_parity] = parity
        c[self.rs.n_parity:] = msg
        return c


_FAST_RS_CACHE = {}


def _get_fast_rs(rs: ShortenedRSCode) -> _FastRSParity:
    key = (rs.m, rs.n, rs.k)
    if key not in _FAST_RS_CACHE:
        _FAST_RS_CACHE[key] = _FastRSParity(rs)
    return _FAST_RS_CACHE[key]


class InterleavedCascadeCodec:
    """RS+BCH cascade with matrix interleaving, x RS codewords per frame."""

    def __init__(self, mode: str = "none", bch_decoder: str = "direct",
                 x_rs: int = PAPER_X_RS):
        assert bch_decoder in ("conv", "direct")
        self.mode = mode
        self.bch_decoder_name = bch_decoder
        self.x_rs = x_rs

        self.m_rs = PAPER_M_RS
        self.n_rs = PAPER_N_RS
        self.k_rs = PAPER_K_RS
        self.rs = ShortenedRSCode(m=self.m_rs, n_short=self.n_rs, k_short=self.k_rs)
        self.rs_enc = _get_fast_rs(self.rs)
        self.bch = BinaryBCH(m=PAPER_M_BCH, t=PAPER_T_BCH,
                             n_short=PAPER_N_BCH, k_short=PAPER_K_BCH)

        self.rs_bits = x_rs * self.n_rs * self.m_rs
        assert self.rs_bits % self.bch.k == 0, "frame must tile BCH blocks exactly"
        self.n_bch_blocks = self.rs_bits // self.bch.k
        self.n_coded_bits = self.n_bch_blocks * self.bch.n

        self.n_info_symbols = x_rs * self.k_rs
        self.n_info_bits = self.n_info_symbols * self.m_rs
        self.effective_rate = self.n_info_bits / self.n_coded_bits

        self.il = MatrixInterleaver(self.n_bch_blocks, self.bch.n, self.m_rs, mode)

    # ------------------------------------------------------------------
    def encode(self, msg_symbols: np.ndarray) -> np.ndarray:
        """msg_symbols: x_rs * k_rs symbols (concatenated per RS codeword)."""
        assert msg_symbols.size == self.n_info_symbols
        codewords = []
        for c in range(self.x_rs):
            m_c = msg_symbols[c * self.k_rs:(c + 1) * self.k_rs]
            codewords.append(self.rs_enc.encode_systematic(m_c))   # n_rs symbols

        # Step-1: proposed schemes spread symbols across RS codewords.
        if self.il.uses_step1():
            sym_stream = rs_symbol_interleave(codewords)        # x*n_rs symbols
        else:
            sym_stream = np.concatenate(codewords)

        # serialise symbols -> bits (LSB first, matching v3), vectorised
        bit_idx = np.arange(self.m_rs)
        rs_bits = ((sym_stream[:, None] >> bit_idx) & 1).astype(np.int8).reshape(-1)

        # BCH inner encode each block
        coded = np.zeros(self.n_coded_bits, dtype=np.int8)
        for b in range(self.n_bch_blocks):
            block = rs_bits[b * self.bch.k:(b + 1) * self.bch.k]
            coded[b * self.bch.n:(b + 1) * self.bch.n] = self.bch.encode(block)

        # Steps 2-3: transmit-order interleave
        return self.il.interleave(coded)

    def decode(self, rx_bits: np.ndarray, counters: OpCounters = None):
        if counters is None:
            counters = OpCounters()

        # invert Steps 2-3
        coded = self.il.deinterleave(rx_bits)

        # BCH inner decode
        rec_bits = np.zeros(self.n_bch_blocks * self.bch.k, dtype=np.int8)
        n_bch_ok = 0
        for b in range(self.n_bch_blocks):
            block = coded[b * self.bch.n:(b + 1) * self.bch.n]
            if self.bch_decoder_name == "conv":
                cw_hat, ok = self.bch.decode_conventional(block, counters)
            else:
                cw_hat, ok = self.bch.decode_direct(block, counters)
            n_bch_ok += int(ok)
            rec_bits[b * self.bch.k:(b + 1) * self.bch.k] = \
                self.bch.extract_message(cw_hat)

        # bits -> symbols (LSB first), vectorised
        n_sym = self.x_rs * self.n_rs
        bit_idx = np.arange(self.m_rs)
        bits_mat = rec_bits.reshape(n_sym, self.m_rs).astype(np.int64)
        sym_stream = (bits_mat << bit_idx).sum(axis=1).astype(np.int64)

        # invert Step-1
        if self.il.uses_step1():
            codewords = rs_symbol_deinterleave(sym_stream, self.x_rs)
        else:
            codewords = [sym_stream[c * self.n_rs:(c + 1) * self.n_rs]
                         for c in range(self.x_rs)]

        # RS outer decode each codeword
        msg_hat = np.zeros(self.n_info_symbols, dtype=np.int64)
        n_rs_ok = 0
        for c in range(self.x_rs):
            c_dec, ok = self.rs.bm_decode(codewords[c], counters)
            n_rs_ok += int(ok)
            if ok:
                mh = self.rs.extract_message(c_dec)
            else:
                mh = codewords[c][self.n_rs - self.k_rs:]
            msg_hat[c * self.k_rs:(c + 1) * self.k_rs] = mh

        return msg_hat, (n_rs_ok == self.x_rs), {
            "n_bch_ok": n_bch_ok, "n_rs_ok": n_rs_ok, "counters": counters}


# =====================================================================
# Channel dispatch (AWGN or burst)
# =====================================================================
def run_channel(bits, ebn0_db, rate, p, rng):
    """p == 0 -> memoryless AWGN; p > 0 -> DFE error-propagation burst channel."""
    if p <= 0.0:
        return run_pam4_channel_hard(bits, ebn0_db, rate, rng)
    return run_pam4_channel_burst(bits, ebn0_db, rate, p, rng)
