"""PAM4 modulator + AWGN channel + per-bit LLR calculator.

Constellation (Gray-encoded, 每符号 2 bits: b1 b0):
    电平 | b1 b0
    -3   | 0 0
    -1   | 0 1
    +1   | 1 1
    +3   | 1 0

Note: b1 is the MSB (upper half), b0 is the LSB (inner).

平均符号能量 E_s = (9+1+1+9)/4 = 5
每比特能量 E_b = E_s / 2 = 2.5 (PAM4 每符号 2 bits)

参考:
- Cioffi, EE379, PAM decoder handout.
- Xiong, Digital Modulation Techniques, Chap 3.
"""
from __future__ import annotations

import numpy as np


# --- Constellation constants -------------------------------------------
# levels[i] = symbol value for bit pair (b1 b0) encoded as int i=2*b1+b0
# But Gray encoding: (b1 b0) -> (00, 01, 11, 10) = levels (-3, -1, +1, +3)
PAM4_LEVELS = np.array([-3.0, -1.0, +3.0, +1.0])
# ↑ index 0 = 00 = -3, index 1 = 01 = -1, index 2 = 10 = +3, index 3 = 11 = +1

# Alternative easier lookup: for bit pair (b1, b0), symbol is:
# (b1, b0) = (0,0) -> -3
# (b1, b0) = (0,1) -> -1
# (b1, b0) = (1,1) -> +1
# (b1, b0) = (1,0) -> +3
BITS_TO_LEVEL = {
    (0, 0): -3.0,
    (0, 1): -1.0,
    (1, 1): +1.0,
    (1, 0): +3.0,
}

# All 4 symbols
ALL_LEVELS = np.array([-3.0, -1.0, +1.0, +3.0])
# For each level, its (b1, b0)
LEVEL_TO_BITS = {
    -3.0: (0, 0),
    -1.0: (0, 1),
    +1.0: (1, 1),
    +3.0: (1, 0),
}

E_S_AVG = 5.0        # average symbol energy
E_B_AVG = E_S_AVG / 2.0  # average bit energy


def bits_to_pam4(bits: np.ndarray) -> np.ndarray:
    """Map a binary vector (length 2N) to N PAM4 symbols.
    bits[2*i]   is b1 (MSB), bits[2*i+1] is b0 (LSB).
    """
    assert bits.size % 2 == 0
    n_syms = bits.size // 2
    out = np.zeros(n_syms)
    for i in range(n_syms):
        b1 = int(bits[2 * i])
        b0 = int(bits[2 * i + 1])
        out[i] = BITS_TO_LEVEL[(b1, b0)]
    return out


def pam4_to_bits_hard(y: np.ndarray) -> np.ndarray:
    """Hard-decision demodulator: nearest level → bit pair."""
    out = np.zeros(2 * y.size, dtype=np.int8)
    for i in range(y.size):
        # nearest of {-3, -1, +1, +3}
        level = ALL_LEVELS[np.argmin(np.abs(ALL_LEVELS - y[i]))]
        b1, b0 = LEVEL_TO_BITS[level]
        out[2 * i] = b1
        out[2 * i + 1] = b0
    return out


def sigma_from_ebn0_pam4(ebn0_db: float, rate: float) -> float:
    """Given Eb/N0 in dB and code rate, return sigma for AWGN.

    E_b (raw, per info bit) = E_b_avg (PAM4 avg) * rate  ...
    Actually the convention: SNR is measured at the receiver where each
    channel bit has energy E_b_avg = 2.5 * rate (accounting for coding).

    Simpler: define
        sigma^2 = E_s / (2 * SNR_per_symbol)
    where SNR_per_symbol = (Eb/N0) * (2 * rate)  (since 2 bits/symbol, rate
    is info-to-coded ratio).

    We follow the common convention: `Eb/N0` refers to information bit energy.
    Each transmitted PAM4 symbol carries 2 coded bits = 2/rate info bits.
    So E_b(info) = E_s / (2/rate) = E_s * rate / 2.
    N0 = 2*sigma^2 (double-sided PSD).

    Solve: sigma^2 = E_s * rate / (2 * 2 * Eb/N0) = E_s * rate / (4 * Eb/N0)
    """
    ebn0_lin = 10.0 ** (ebn0_db / 10.0)
    sigma2 = E_S_AVG * rate / (4.0 * ebn0_lin)
    return float(np.sqrt(sigma2))


def awgn_channel(x: np.ndarray, sigma: float, rng: np.random.Generator) -> np.ndarray:
    return x + sigma * rng.standard_normal(x.shape)


def pam4_bit_llr(y: np.ndarray, sigma: float, use_max_log: bool = False) -> np.ndarray:
    """Compute per-bit LLR for PAM4 received symbols.

    Returns a flat array of length 2N where output[2*i]   = LLR(b1|y_i)
    and output[2*i+1] = LLR(b0|y_i).

    LLR convention: LLR = log P(b=0|y) / P(b=1|y).
    Positive LLR => bit 0 is more likely.

    Uses full log-sum-exp by default; max-log-approx when use_max_log=True.
    """
    inv_2s2 = 1.0 / (2.0 * sigma * sigma)
    # For each candidate level s, compute exponent -||y-s||^2 / (2 sigma^2)
    n = y.size
    out = np.zeros(2 * n)
    for i in range(n):
        # b1 = 0: levels -3, -1 (i.e., 00 and 01)
        # b1 = 1: levels +1, +3 (i.e., 11 and 10)
        d_m3 = -inv_2s2 * (y[i] - (-3.0)) ** 2
        d_m1 = -inv_2s2 * (y[i] - (-1.0)) ** 2
        d_p1 = -inv_2s2 * (y[i] - (+1.0)) ** 2
        d_p3 = -inv_2s2 * (y[i] - (+3.0)) ** 2

        # LLR(b1|y) = log( exp(d_m3)+exp(d_m1) ) / ( exp(d_p1)+exp(d_p3) )
        # LLR(b0|y) = log( exp(d_m3)+exp(d_p3) ) / ( exp(d_m1)+exp(d_p1) )
        if use_max_log:
            llr_b1 = max(d_m3, d_m1) - max(d_p1, d_p3)
            llr_b0 = max(d_m3, d_p3) - max(d_m1, d_p1)
        else:
            llr_b1 = _logsumexp2(d_m3, d_m1) - _logsumexp2(d_p1, d_p3)
            llr_b0 = _logsumexp2(d_m3, d_p3) - _logsumexp2(d_m1, d_p1)
        out[2 * i] = llr_b1
        out[2 * i + 1] = llr_b0
    return out


def _logsumexp2(a: float, b: float) -> float:
    """log(exp(a) + exp(b)) - numerically stable."""
    m = a if a > b else b
    return m + np.log1p(np.exp(-abs(a - b)))
