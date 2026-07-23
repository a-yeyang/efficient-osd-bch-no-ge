"""Monte Carlo BER-SNR simulation runner for the cascade codec.

Runs pure RS BM, pure RS LCC-BR, cascade scheme A, cascade scheme B
across an SNR sweep and reports:
- Symbol error rate (SER) on message symbols
- Bit error rate (BER)
- Frame error rate (FER)
- Per-frame decoding ops (F_{2^m}) and clock cycles
"""
from __future__ import annotations

import numpy as np
import time
from dataclasses import dataclass, field, asdict
from typing import Callable
import json

from .cascade import CascadeConfig, CascadedCodec, PureRSCodec, run_channel
from .lagrange_cache import LagrangeCache, cascade_scheme_b_decode
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
    avg_latency_us: list = field(default_factory=list)


def _bits_of_symbols(syms: np.ndarray, m: int) -> np.ndarray:
    """Convert symbol array (n values, each < 2^m) to bit array (n*m bits)."""
    n = syms.size
    out = np.zeros(n * m, dtype=np.int8)
    for i in range(n):
        for b in range(m):
            out[i * m + b] = (int(syms[i]) >> b) & 1
    return out


def run_bench(
    method_name: str,
    codec_pack,  # tuple: (encoder_fn, decoder_fn, effective_rate)
    ebn0_list,
    n_info_symbols: int,
    m_bits_per_symbol: int,
    seed: int = 0,
    min_frame_errors: int = 30,
    max_frames: int = 5000,
    verbose: bool = True,
):
    """Run MC sim for a given codec across SNR."""
    encoder, decoder, rate = codec_pack
    zero_msg = np.zeros(n_info_symbols, dtype=np.int64)

    res = SimResult()
    for ebn0 in ebn0_list:
        rng = np.random.default_rng(seed + int(round(ebn0 * 100)))
        n_frames = 0
        n_frame_errors = 0
        n_bit_errors = 0
        n_symbol_errors = 0
        sum_f2m = 0.0
        sum_lat = 0.0
        t_start = time.perf_counter()

        while n_frames < max_frames:
            # Encode all-zero for symmetry (BPSK/PAM4 Gray symmetric).
            # Actually PAM4 Gray isn't fully symmetric. Use random msg instead.
            msg = rng.integers(0, 1 << m_bits_per_symbol, n_info_symbols)
            coded = encoder(msg)
            llr = run_channel(coded, ebn0, rate, rng)

            t_dec_start = time.perf_counter()
            counters = OpCounters()
            msg_hat, res_info = decoder(llr, counters)
            t_dec = (time.perf_counter() - t_dec_start) * 1e6

            n_frames += 1
            # Symbol errors
            n_sym_err_this = int(np.sum(msg_hat != msg))
            if n_sym_err_this > 0:
                n_frame_errors += 1
            n_symbol_errors += n_sym_err_this
            # Bit errors
            msg_bits = _bits_of_symbols(msg, m_bits_per_symbol)
            msg_hat_bits = _bits_of_symbols(msg_hat, m_bits_per_symbol)
            n_bit_errors += int(np.sum(msg_bits != msg_hat_bits))
            sum_f2m += counters.f2m
            sum_lat += t_dec
            if n_frame_errors >= min_frame_errors and n_frames >= 100:
                break

        elapsed = time.perf_counter() - t_start
        ser = n_symbol_errors / max(1, n_frames * n_info_symbols)
        ber = n_bit_errors / max(1, n_frames * n_info_symbols * m_bits_per_symbol)
        fer = n_frame_errors / max(1, n_frames)
        avg_f2m = sum_f2m / max(1, n_frames)
        avg_lat = sum_lat / max(1, n_frames)

        if verbose:
            print(f"  {method_name} @ {ebn0:.2f} dB: FER={fer:.3e}, "
                  f"BER={ber:.3e}, avg_f2m={avg_f2m:.0f}, "
                  f"lat={avg_lat:.0f}μs, "
                  f"{n_frames} frames, {elapsed:.1f}s")

        res.ebn0_db.append(ebn0)
        res.ser.append(ser)
        res.ber.append(ber)
        res.fer.append(fer)
        res.n_frames.append(n_frames)
        res.n_frame_errors.append(n_frame_errors)
        res.avg_f2m_ops.append(avg_f2m)
        res.avg_latency_us.append(avg_lat)

        if fer < 1e-6 and n_frame_errors < 3:
            break

    return res


def get_codec_pack(cfg: CascadeConfig, method: str):
    """Return (encoder_fn, decoder_fn, effective_rate) for a given method."""
    if method in ("pure_rs_bm", "pure_rs_lccbr"):
        codec = PureRSCodec(cfg)
        encoder = codec.encode
        rs_method = "hard" if method == "pure_rs_bm" else "soft"

        def decoder(llr, counters):
            msg_hat, res = codec.decode(llr, method=rs_method, counters=counters)
            return msg_hat, res
        return encoder, decoder, codec.effective_rate

    else:  # cascade
        codec = CascadedCodec(cfg)
        encoder = codec.encode

        if method == "scheme_a":
            def decoder(llr, counters):
                return codec.decode(llr, method="scheme_a", counters=counters)
        elif method == "scheme_b":
            def decoder(llr, counters):
                return cascade_scheme_b_decode(codec, llr, counters=counters)
        else:
            raise ValueError(method)
        return encoder, decoder, codec.effective_rate


def result_to_dict(res: SimResult) -> dict:
    return asdict(res)
