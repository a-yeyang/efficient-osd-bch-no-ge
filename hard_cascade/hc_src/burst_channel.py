"""Burst-error channel for the interleaving study (advisor's 4th task, 2026-07-27).

Reproduces the error-propagation channel of Liu/Song/Wang, "A Novel Interleaving
Scheme for Concatenated Codes on Burst-Error Channel," APCC 2022.

The paper models burst errors with a one-tap Decision-Feedback Equalizer (DFE):
the PAM4 symbols are convolved with a one-tap ISI channel, AWGN is added, then a
DFE recovers the symbols in the receiver.  A wrong decision feeds back and makes
the *next* symbol more likely to be wrong, producing bursts.  The paper
parameterises this by

    p = P(next signal wrong | current signal wrong)          (error propagation)

and sweeps p in {0.3, 0.5, 0.75}.  p = 0 is the memoryless AWGN channel.

We reproduce the *statistics* the paper cares about (burst runs at the PAM4
*symbol* level with geometric length ~1/(1-p)) with a faithful and reproducible
model:

  * Each PAM4 symbol first errs independently from AWGN (the "seed" error), with
    the same hard-decision demapper used everywhere else in the project.
  * Whenever a symbol is in error (seed OR propagated), the next symbol is forced
    to be in error with probability p (DFE feedback of a wrong decision).  A
    forced symbol error flips the symbol to its nearest constellation neighbour,
    which changes 1 of its 2 Gray-coded bits (the realistic DFE behaviour — a
    slightly-off decision slips to an adjacent level, a single Gray-bit flip).

This keeps p = 0 bit-identical in distribution to the existing AWGN path, while
p > 0 injects geometric symbol-error bursts exactly as the paper's DFE does.
"""
from __future__ import annotations

import numpy as np

from .upstream import (
    bits_to_pam4, pam4_to_bits_hard, sigma_from_ebn0_pam4, awgn_channel,
)

# PAM4 Gray levels in ascending order and their (b1,b0) mapping (see pam4.py).
_ASC_LEVELS = np.array([-3.0, -1.0, +1.0, +3.0])
_LEVEL_BITS = {-3.0: (0, 0), -1.0: (0, 1), +1.0: (1, 1), +3.0: (1, 0)}
# ascending-index -> (b1, b0), for vectorised symbol->bit mapping
_LEVEL_BITS_ARR = np.array([[0, 0], [0, 1], [1, 1], [1, 0]], dtype=np.int8)


def _nearest_neighbour_level(level_idx: int, rng: np.random.Generator) -> int:
    """Index of an adjacent constellation level (a DFE slip to a neighbour).

    Interior levels slip up or down with equal probability; edge levels slip
    inward.  Adjacent Gray levels differ in exactly one bit, so this is a
    single-bit symbol error — the realistic 'slightly wrong decision' burst.
    """
    if level_idx == 0:
        return 1
    if level_idx == len(_ASC_LEVELS) - 1:
        return len(_ASC_LEVELS) - 2
    return level_idx + (1 if rng.random() < 0.5 else -1)


def run_pam4_channel_burst(bits: np.ndarray, ebn0_db: float, rate: float,
                           p: float, rng: np.random.Generator) -> np.ndarray:
    """Hard-decision PAM4 over an error-propagation (burst) channel.

    Parameters
    ----------
    bits : transmitted coded bits (0/1), length even after internal padding.
    ebn0_db, rate : as for the AWGN path (seed error rate).
    p : error-propagation probability P(next sym wrong | current sym wrong).
        p = 0 recovers the memoryless AWGN channel exactly.
    rng : numpy Generator.

    Returns the hard-decision received bits (same length as `bits`).
    """
    assert 0.0 <= p < 1.0
    n_orig = bits.size
    if n_orig % 2 != 0:
        bits_padded = np.concatenate([bits, np.zeros(1, dtype=np.int8)])
    else:
        bits_padded = bits

    x = bits_to_pam4(bits_padded)
    sigma = sigma_from_ebn0_pam4(ebn0_db, rate)
    y = awgn_channel(x, sigma, rng)

    # Step 1: memoryless hard decision (seed errors from AWGN).
    # nearest of {-3,-1,+1,+3} via decision thresholds at {-2,0,+2} (vectorised).
    tx_idx = ((x + 3.0) / 2.0).round().astype(np.int64)          # exact (on-grid)
    rx_idx = np.digitize(y, [-2.0, 0.0, 2.0]).astype(np.int64)

    if p > 0.0:
        # Step 2: DFE error propagation. Walk symbols left→right; if the current
        # symbol is wrong, force the next wrong with probability p.
        n_sym = rx_idx.size
        u = rng.random(n_sym)
        for i in range(n_sym - 1):
            cur_wrong = rx_idx[i] != tx_idx[i]
            if cur_wrong and rx_idx[i + 1] == tx_idx[i + 1] and u[i + 1] < p:
                # force next symbol into a single-Gray-bit-flip neighbour error
                rx_idx[i + 1] = _nearest_neighbour_level(int(tx_idx[i + 1]), rng)

    # Symbols → bits via the Gray map (vectorised).
    out = _LEVEL_BITS_ARR[rx_idx].reshape(-1)
    return out[:n_orig].astype(np.int8)


def measure_burst_stats(ebn0_db: float, rate: float, p: float,
                        n_sym: int = 200_000, seed: int = 7):
    """Empirical symbol-error rate and mean burst-run length, for validation."""
    rng = np.random.default_rng(seed)
    bits = rng.integers(0, 2, 2 * n_sym).astype(np.int8)
    x = bits_to_pam4(bits)
    sigma = sigma_from_ebn0_pam4(ebn0_db, rate)
    y = awgn_channel(x, sigma, rng)
    tx_idx = np.searchsorted(_ASC_LEVELS, x)
    rx_idx = np.array([int(np.argmin(np.abs(_ASC_LEVELS - yi))) for yi in y])
    u = rng.random(n_sym)
    for i in range(n_sym - 1):
        if rx_idx[i] != tx_idx[i] and rx_idx[i + 1] == tx_idx[i + 1] and u[i + 1] < p:
            rx_idx[i + 1] = _nearest_neighbour_level(int(tx_idx[i + 1]), rng)
    wrong = rx_idx != tx_idx
    ser = wrong.mean()
    # mean run length of consecutive wrong symbols
    runs = []
    run = 0
    for w in wrong:
        if w:
            run += 1
        elif run:
            runs.append(run)
            run = 0
    if run:
        runs.append(run)
    mean_run = float(np.mean(runs)) if runs else 0.0
    return {"ser": float(ser), "mean_run": mean_run,
            "n_bursts": len(runs), "expected_run_if_geom": 1.0 / (1.0 - p) if p < 1 else np.inf}
