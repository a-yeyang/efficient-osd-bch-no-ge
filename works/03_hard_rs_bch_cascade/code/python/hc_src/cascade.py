"""Hard-decision RS+BCH cascade codec.

Encoding chain:
    msg (k_RS symbols) → RS encode → n_RS symbols → serialize to (m * n_RS) bits
    → partition into k_BCH-sized chunks → BCH encode each → concatenate
    → PAM4 modulate → AWGN

Decoding chain:
    y (PAM4 symbols) → hard-decision demod → bit stream
    → partition into n_BCH-sized chunks → BCH decode (Conv or Direct) → extract k_BCH msg
    → re-assemble to (m * n_RS) bit stream → parse into n_RS symbols
    → RS BM decode → extract k_RS message

Because BCH encodes bits and RS encodes symbols, we need careful bit/symbol packing.
"""
from __future__ import annotations

import numpy as np
from dataclasses import dataclass

from .upstream import (
    OpCounters, RSCode,
    bits_to_pam4, pam4_to_bits_hard, sigma_from_ebn0_pam4, awgn_channel,
)
from .bch_t2 import BCHt2Code
from .latency_model import LatencyModel


@dataclass
class HardCascadeConfig:
    m: int
    k_rs: int          # RS message length in symbols
    # BCH is BCHt2Code(m), so t_bch = 2 fixed

    @property
    def n_rs(self) -> int:
        return (1 << self.m) - 1

    @property
    def t_rs(self) -> int:
        return (self.n_rs - self.k_rs) // 2

    def describe(self) -> str:
        n = self.n_rs
        return f"RS({n},{self.k_rs}, t={self.t_rs}) + BCH({n}, t=2)"


class HardCascadedCodec:
    """RS + BCH t=2 hard-decision cascade."""

    def __init__(self, cfg: HardCascadeConfig, bch_decoder: str = "direct"):
        """bch_decoder: 'conv' (BM+Chien) or 'direct' (LUT root finding)."""
        self.cfg = cfg
        self.rs = RSCode(m=cfg.m, k=cfg.k_rs)
        self.bch = BCHt2Code(m=cfg.m)
        assert bch_decoder in ("conv", "direct")
        self.bch_decoder_name = bch_decoder

        self.rs_bits = cfg.m * cfg.n_rs  # bits after RS encoding
        # BCH encodes k_BCH bits at a time -> n_BCH bits.
        # We need to pad rs_bits to a multiple of k_BCH.
        if self.rs_bits % self.bch.k != 0:
            self.n_pad_bits = self.bch.k - (self.rs_bits % self.bch.k)
        else:
            self.n_pad_bits = 0
        self.n_bch_blocks = (self.rs_bits + self.n_pad_bits) // self.bch.k
        self.n_coded_bits = self.n_bch_blocks * self.bch.n

        # Rate accounting for SNR translation
        self.n_info_bits = cfg.k_rs * cfg.m
        self.effective_rate = self.n_info_bits / self.n_coded_bits

    # -----------------------------------------------------
    def encode(self, msg_symbols: np.ndarray) -> np.ndarray:
        """msg (k_RS symbols) → cascade bits."""
        cfg = self.cfg
        assert msg_symbols.size == cfg.k_rs
        c_rs = self.rs.encode_systematic(msg_symbols)  # n_rs symbols
        # Serialize to bits (LSB first)
        rs_bits = np.zeros(self.rs_bits, dtype=np.int8)
        for i, s in enumerate(c_rs):
            for b in range(cfg.m):
                rs_bits[i * cfg.m + b] = (int(s) >> b) & 1
        # Pad
        if self.n_pad_bits > 0:
            rs_bits = np.concatenate([rs_bits, np.zeros(self.n_pad_bits, dtype=np.int8)])
        # BCH encode each block
        out = np.zeros(self.n_coded_bits, dtype=np.int8)
        for b in range(self.n_bch_blocks):
            block = rs_bits[b * self.bch.k:(b + 1) * self.bch.k]
            out[b * self.bch.n:(b + 1) * self.bch.n] = self.bch.encode(block)
        return out

    def decode(self, hard_bits: np.ndarray, counters: OpCounters = None):
        """Hard-decision decode. Returns (msg_hat, ok, stats)."""
        cfg = self.cfg
        if counters is None:
            counters = OpCounters()

        # 1. BCH inner decode each block
        recovered_bits = np.zeros(self.n_bch_blocks * self.bch.k, dtype=np.int8)
        n_bch_ok = 0
        for b in range(self.n_bch_blocks):
            block = hard_bits[b * self.bch.n:(b + 1) * self.bch.n]
            if self.bch_decoder_name == "conv":
                cw_hat, ok = self.bch.decode_conventional(block, counters)
            else:
                cw_hat, ok = self.bch.decode_direct(block, counters)
            if ok:
                n_bch_ok += 1
            # Extract message: since our BCH is systematic with parity first,
            # message = cw_hat[n-k:] (last k bits)
            recovered_bits[b * self.bch.k:(b + 1) * self.bch.k] = self.bch.extract_message(cw_hat)

        # 2. Reassemble into RS symbols
        rs_bit_stream = recovered_bits[:self.rs_bits]
        r_rs = np.zeros(cfg.n_rs, dtype=np.int64)
        for i in range(cfg.n_rs):
            s = 0
            for b in range(cfg.m):
                s |= (int(rs_bit_stream[i * cfg.m + b]) & 1) << b
            r_rs[i] = s

        # 3. RS BM decode
        c_dec, ok_rs = self.rs.bm_decode(r_rs, counters)
        if ok_rs:
            msg_hat = self.rs.extract_message(c_dec)
        else:
            msg_hat = r_rs[cfg.n_rs - cfg.k_rs:]

        return msg_hat, ok_rs, {
            "n_bch_ok": n_bch_ok,
            "counters": counters,
        }

    # -----------------------------------------------------
    def latency_cycles(self, lagrange_mode: str = "none") -> int:
        """Total cascade clock cycles.

        lagrange_mode:
            'none'      no sharing.
            'v1'        share alpha table + syndrome unit  (save 1 cyc)
            'v2'        v1 + share GF multiplier array     (save 2 cyc)
            True        (backward-compat) same as 'v1'
            False       (backward-compat) same as 'none'
        """
        # BCH cycles depend on decoder AND (for Direct) on field size m
        if self.bch_decoder_name == "conv":
            bch_cyc = LatencyModel.bch_conv_cycles(2, m=self.cfg.m)
        else:
            bch_cyc = LatencyModel.bch_direct_cycles(2, m=self.cfg.m)

        # Normalize legacy boolean
        if lagrange_mode is True:
            lagrange_mode = "v1"
        elif lagrange_mode is False:
            lagrange_mode = "none"

        if lagrange_mode == "none":
            return LatencyModel.cascade_serial(bch_cyc, self.cfg.t_rs)
        elif lagrange_mode == "v1":
            return LatencyModel.cascade_lagrange_v1(bch_cyc, self.cfg.t_rs)
        elif lagrange_mode == "v2":
            return LatencyModel.cascade_lagrange_v2(bch_cyc, self.cfg.t_rs)
        else:
            raise ValueError(f"unknown lagrange_mode {lagrange_mode}")


class PureRSHardCodec:
    """Pure RS baseline (no BCH inner): hard-decision decoding only."""

    def __init__(self, cfg: HardCascadeConfig):
        self.cfg = cfg
        self.rs = RSCode(m=cfg.m, k=cfg.k_rs)
        self.rs_bits = cfg.m * cfg.n_rs
        self.n_coded_bits = self.rs_bits
        self.n_info_bits = cfg.k_rs * cfg.m
        self.effective_rate = self.n_info_bits / self.n_coded_bits

    def encode(self, msg_symbols: np.ndarray) -> np.ndarray:
        cfg = self.cfg
        c_rs = self.rs.encode_systematic(msg_symbols)
        bits = np.zeros(cfg.m * cfg.n_rs, dtype=np.int8)
        for i, s in enumerate(c_rs):
            for b in range(cfg.m):
                bits[i * cfg.m + b] = (int(s) >> b) & 1
        return bits

    def decode(self, hard_bits: np.ndarray, counters: OpCounters = None):
        cfg = self.cfg
        if counters is None:
            counters = OpCounters()

        r_rs = np.zeros(cfg.n_rs, dtype=np.int64)
        for i in range(cfg.n_rs):
            s = 0
            for b in range(cfg.m):
                s |= (int(hard_bits[i * cfg.m + b]) & 1) << b
            r_rs[i] = s

        c_dec, ok = self.rs.bm_decode(r_rs, counters)
        if ok:
            msg_hat = self.rs.extract_message(c_dec)
        else:
            msg_hat = r_rs[cfg.n_rs - cfg.k_rs:]
        return msg_hat, ok, {"counters": counters}

    def latency_cycles(self) -> int:
        return LatencyModel.rs_bm_cycles(self.cfg.t_rs)


# =====================================================================
# PAM4 channel wrapper (hard-decision only)
# =====================================================================
def run_pam4_channel_hard(bits: np.ndarray, ebn0_db: float, rate: float,
                          rng: np.random.Generator) -> np.ndarray:
    """Modulate PAM4 → AWGN → hard-decision demod. Returns hard bits."""
    sigma = sigma_from_ebn0_pam4(ebn0_db, rate)
    n_orig = bits.size
    if n_orig % 2 != 0:
        bits_padded = np.concatenate([bits, np.zeros(1, dtype=np.int8)])
    else:
        bits_padded = bits
    x = bits_to_pam4(bits_padded)
    y = awgn_channel(x, sigma, rng)
    hard = pam4_to_bits_hard(y)
    return hard[:n_orig]
