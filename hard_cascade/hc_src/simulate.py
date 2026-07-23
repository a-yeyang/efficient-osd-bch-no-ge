"""Monte Carlo simulation driver for hard-decision cascade.

Runs BER/FER + latency for:
  - Pure RS-BM (baseline)
  - Cascade Conv (BCH conventional + RS-BM)
  - Cascade Direct (BCH direct + RS-BM)
  - Cascade Conv + Lagrange share
  - Cascade Direct + Lagrange share
"""
from __future__ import annotations

import time
import numpy as np
from dataclasses import dataclass, field, asdict
from typing import Callable

from .cascade import (
    HardCascadeConfig, HardCascadedCodec, PureRSHardCodec, run_pam4_channel_hard,
)
from .upstream import OpCounters


@dataclass
class SimResult:
    ebn0_db: list = field(default_factory=list)
    ser: list = field(default_factory=list)
    ber: list = field(default_factory=list)
    fer: list = field(default_factory=list)
    n_frames: list = field(default_factory=list)
    n_frame_errors: list = field(default_factory=list)
    avg_f2m_ops: list = field(default_factory=list)


def _bits_of_symbols(syms: np.ndarray, m: int) -> np.ndarray:
    n = syms.size
    out = np.zeros(n * m, dtype=np.int8)
    for i in range(n):
        for b in range(m):
            out[i * m + b] = (int(syms[i]) >> b) & 1
    return out


def run_bench(
    method_name: str,
    codec,
    ebn0_list,
    n_info_symbols: int,
    m_bits_per_symbol: int,
    seed: int = 0,
    min_frame_errors: int = 30,
    max_frames: int = 2000,
    verbose: bool = True,
) -> SimResult:
    """Monte Carlo over SNR points."""
    res = SimResult()
    for ebn0 in ebn0_list:
        rng = np.random.default_rng(seed + int(round(ebn0 * 100)))
        n_frames = 0
        n_frame_errors = 0
        n_bit_errors = 0
        n_sym_errors = 0
        sum_f2m = 0.0
        t_start = time.perf_counter()
        while n_frames < max_frames:
            msg = rng.integers(0, 1 << m_bits_per_symbol, n_info_symbols)
            coded = codec.encode(msg)
            hard = run_pam4_channel_hard(coded, ebn0, codec.effective_rate, rng)
            counters = OpCounters()
            msg_hat, _ok, _ = codec.decode(hard, counters)

            n_frames += 1
            n_this = int(np.sum(msg_hat != msg))
            if n_this > 0:
                n_frame_errors += 1
            n_sym_errors += n_this
            msg_bits = _bits_of_symbols(msg, m_bits_per_symbol)
            msg_hat_bits = _bits_of_symbols(msg_hat, m_bits_per_symbol)
            n_bit_errors += int(np.sum(msg_bits != msg_hat_bits))
            sum_f2m += counters.f2m
            if n_frame_errors >= min_frame_errors and n_frames >= 100:
                break

        ser = n_sym_errors / max(1, n_frames * n_info_symbols)
        ber = n_bit_errors / max(1, n_frames * n_info_symbols * m_bits_per_symbol)
        fer = n_frame_errors / max(1, n_frames)
        avg_f2m = sum_f2m / max(1, n_frames)
        dt = time.perf_counter() - t_start

        if verbose:
            print(f"  {method_name} @ {ebn0:.2f} dB: FER={fer:.3e}, "
                  f"BER={ber:.3e}, avg_f2m={avg_f2m:.0f}, "
                  f"{n_frames} frames, {dt:.1f}s")

        res.ebn0_db.append(ebn0)
        res.ser.append(ser)
        res.ber.append(ber)
        res.fer.append(fer)
        res.n_frames.append(n_frames)
        res.n_frame_errors.append(n_frame_errors)
        res.avg_f2m_ops.append(avg_f2m)

        if fer < 1e-6 and n_frame_errors < 2:
            break

    return res


def result_to_dict(res: SimResult) -> dict:
    return asdict(res)
