"""v3 cascade codec: decoupled RS field and BCH field, shortened codes.

Supports the advisor's third-task configs:

  Config 1: RS(544,514, t=15)/GF(2^10)  +  BCH(144,136, t=1)/GF(2^8)
  Config 2: RS(255,239, t=8)/GF(2^8)    +  BCH(255,239, t=2)/GF(2^8)

Unlike the v1 HardCascadedCodec (which tied RS and BCH to one field m), here
the outer RS field (m_rs) and inner BCH field (m_bch) are independent.

Chain:
  msg (k_rs symbols /GF(2^m_rs)) -> RS encode -> n_rs symbols
    -> serialise to m_rs * n_rs bits (LSB first)
    -> partition into k_bch-bit blocks (zero-pad the tail) -> BCH encode each
    -> concatenate -> PAM4 -> AWGN
"""
from __future__ import annotations

import numpy as np
from dataclasses import dataclass

from .upstream import (
    OpCounters,
    bits_to_pam4, pam4_to_bits_hard, sigma_from_ebn0_pam4, awgn_channel,
)
from .shortened_codes import ShortenedRSCode, BinaryBCH
from .latency_model import LatencyModel


@dataclass
class CascadeV3Config:
    name: str
    # RS (outer)
    m_rs: int
    n_rs: int
    k_rs: int
    # BCH (inner)
    m_bch: int
    t_bch: int
    n_bch: int
    k_bch: int

    @property
    def t_rs(self) -> int:
        return (self.n_rs - self.k_rs) // 2

    def describe(self) -> str:
        return (f"{self.name}: RS({self.n_rs},{self.k_rs}, t={self.t_rs})"
                f"/GF(2^{self.m_rs}) + BCH({self.n_bch},{self.k_bch}, "
                f"t={self.t_bch})/GF(2^{self.m_bch})")


CONFIGS_V3 = {
    "cfg1_kp4": CascadeV3Config(
        name="Config1 (KP4)",
        m_rs=10, n_rs=544, k_rs=514,
        m_bch=8, t_bch=1, n_bch=144, k_bch=136,
    ),
    "cfg2_255": CascadeV3Config(
        name="Config2 (255)",
        m_rs=8, n_rs=255, k_rs=239,
        m_bch=8, t_bch=2, n_bch=255, k_bch=239,
    ),
}


class CascadeV3Codec:
    """RS (outer, GF(2^m_rs)) + BCH (inner, GF(2^m_bch)) hard-decision cascade."""

    def __init__(self, cfg: CascadeV3Config, bch_decoder: str = "direct"):
        assert bch_decoder in ("conv", "direct")
        self.cfg = cfg
        self.bch_decoder_name = bch_decoder

        self.rs = ShortenedRSCode(m=cfg.m_rs, n_short=cfg.n_rs, k_short=cfg.k_rs)
        self.bch = BinaryBCH(m=cfg.m_bch, t=cfg.t_bch,
                             n_short=cfg.n_bch, k_short=cfg.k_bch)

        self.rs_bits = cfg.m_rs * cfg.n_rs
        if self.rs_bits % self.bch.k != 0:
            self.n_pad_bits = self.bch.k - (self.rs_bits % self.bch.k)
        else:
            self.n_pad_bits = 0
        self.n_bch_blocks = (self.rs_bits + self.n_pad_bits) // self.bch.k
        self.n_coded_bits = self.n_bch_blocks * self.bch.n

        self.n_info_bits = cfg.k_rs * cfg.m_rs
        self.effective_rate = self.n_info_bits / self.n_coded_bits

    # ------------------------------------------------------------------
    def encode(self, msg_symbols: np.ndarray) -> np.ndarray:
        cfg = self.cfg
        assert msg_symbols.size == cfg.k_rs
        c_rs = self.rs.encode_systematic(msg_symbols)   # n_rs symbols
        rs_bits = np.zeros(self.rs_bits, dtype=np.int8)
        for i, s in enumerate(c_rs):
            si = int(s)
            for b in range(cfg.m_rs):
                rs_bits[i * cfg.m_rs + b] = (si >> b) & 1
        if self.n_pad_bits > 0:
            rs_bits = np.concatenate(
                [rs_bits, np.zeros(self.n_pad_bits, dtype=np.int8)])
        out = np.zeros(self.n_coded_bits, dtype=np.int8)
        for b in range(self.n_bch_blocks):
            block = rs_bits[b * self.bch.k:(b + 1) * self.bch.k]
            out[b * self.bch.n:(b + 1) * self.bch.n] = self.bch.encode(block)
        return out

    def decode(self, hard_bits: np.ndarray, counters: OpCounters = None):
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
            recovered_bits[b * self.bch.k:(b + 1) * self.bch.k] = \
                self.bch.extract_message(cw_hat)

        # 2. Reassemble into RS symbols
        rs_bit_stream = recovered_bits[:self.rs_bits]
        r_rs = np.zeros(cfg.n_rs, dtype=np.int64)
        for i in range(cfg.n_rs):
            s = 0
            for b in range(cfg.m_rs):
                s |= (int(rs_bit_stream[i * cfg.m_rs + b]) & 1) << b
            r_rs[i] = s

        # 3. RS BM decode
        c_dec, ok_rs = self.rs.bm_decode(r_rs, counters)
        if ok_rs:
            msg_hat = self.rs.extract_message(c_dec)
        else:
            msg_hat = r_rs[cfg.n_rs - cfg.k_rs:]

        return msg_hat, ok_rs, {"n_bch_ok": n_bch_ok, "counters": counters}

    # ------------------------------------------------------------------
    def latency_cycles(self, lagrange_mode: str = "none") -> int:
        """Total cascade clock cycles (see latency_model)."""
        if self.bch_decoder_name == "conv":
            bch_cyc = LatencyModel.bch_conv_cycles(self.cfg.t_bch, m=self.cfg.m_bch)
        else:
            bch_cyc = LatencyModel.bch_direct_cycles(self.cfg.t_bch, m=self.cfg.m_bch)
        if lagrange_mode == "none":
            return LatencyModel.cascade_serial(bch_cyc, self.cfg.t_rs)
        elif lagrange_mode == "v1":
            return LatencyModel.cascade_lagrange_v1(bch_cyc, self.cfg.t_rs)
        elif lagrange_mode == "v2":
            return LatencyModel.cascade_lagrange_v2(bch_cyc, self.cfg.t_rs)
        raise ValueError(f"unknown lagrange_mode {lagrange_mode}")


class PureRSV3Codec:
    """Pure RS baseline (no inner BCH), hard-decision only."""

    def __init__(self, cfg: CascadeV3Config):
        self.cfg = cfg
        self.rs = ShortenedRSCode(m=cfg.m_rs, n_short=cfg.n_rs, k_short=cfg.k_rs)
        self.rs_bits = cfg.m_rs * cfg.n_rs
        self.n_coded_bits = self.rs_bits
        self.n_info_bits = cfg.k_rs * cfg.m_rs
        self.effective_rate = self.n_info_bits / self.n_coded_bits

    def encode(self, msg_symbols: np.ndarray) -> np.ndarray:
        cfg = self.cfg
        c_rs = self.rs.encode_systematic(msg_symbols)
        bits = np.zeros(cfg.m_rs * cfg.n_rs, dtype=np.int8)
        for i, s in enumerate(c_rs):
            si = int(s)
            for b in range(cfg.m_rs):
                bits[i * cfg.m_rs + b] = (si >> b) & 1
        return bits

    def decode(self, hard_bits: np.ndarray, counters: OpCounters = None):
        cfg = self.cfg
        if counters is None:
            counters = OpCounters()
        r_rs = np.zeros(cfg.n_rs, dtype=np.int64)
        for i in range(cfg.n_rs):
            s = 0
            for b in range(cfg.m_rs):
                s |= (int(hard_bits[i * cfg.m_rs + b]) & 1) << b
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
# PAM4 channel (hard-decision)
# =====================================================================
def run_pam4_channel_hard(bits: np.ndarray, ebn0_db: float, rate: float,
                          rng: np.random.Generator) -> np.ndarray:
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
