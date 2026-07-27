"""Diagnostic: WHY does no-interleave stay competitive on the burst channel?

Hypothesis: in this exact code geometry one BCH(144,136) block carries only
136/10 = 13.6 RS symbols, which is < the outer RS t=15 budget.  So a short burst
confined to ONE block (the 'none' case) corrupts <=14 symbols in ONE RS codeword
-> RS absorbs it.  Interleaving spreads the burst across many blocks; combined
with the AWGN random-error floor this pushes MORE blocks past t=1, whose scattered
failures corrupt symbols across BOTH RS codewords -> more likely to exceed t=15.

We measure, per frame at a fixed SNR on p=0.75:
  * frame error rate
  * mean # of BCH blocks that fail to decode (n_bch_blocks - n_bch_ok)
  * mean # of RS codewords that fail (x_rs - n_rs_ok)
for modes none / existing / w1bit / w1sym.
"""
from __future__ import annotations
import sys
sys.path.insert(0, ".")
import numpy as np
from hc_src.cascade_v4 import InterleavedCascadeCodec, run_channel

P = 0.75
EBN0 = 10.5
N_FRAMES = 400
MODES = ["none", "existing", "w1sym", "w1bit"]


def bits_of(sym, m):
    return ((sym[:, None] >> np.arange(m)) & 1).astype(np.int8).reshape(-1)


print(f"Diagnostic @ p={P}, Eb/N0={EBN0} dB, {N_FRAMES} frames")
print(f"{'mode':10s} {'FER':>8s} {'bitErr':>9s} {'blkFail':>8s} {'rsFail':>7s}")
for mode in MODES:
    cod = InterleavedCascadeCodec(mode=mode, bch_decoder="direct")
    m = cod.m_rs
    rng = np.random.default_rng(12345 + MODES.index(mode) * 101)
    nfe = blkfail = rsfail = nbe = tot = 0
    for _ in range(N_FRAMES):
        msg = rng.integers(0, 1 << m, cod.n_info_symbols)
        coded = cod.encode(msg)
        rx = run_channel(coded, EBN0, cod.effective_rate, P, rng)
        mh, ok, info = cod.decode(rx)
        e = int((bits_of(msg, m) != bits_of(mh, m)).sum())
        nbe += e
        tot += msg.size * m
        if e > 0:
            nfe += 1
        blkfail += cod.n_bch_blocks - info["n_bch_ok"]
        rsfail += cod.x_rs - info["n_rs_ok"]
    print(f"{mode:10s} {nfe/N_FRAMES:8.4f} {nbe/tot:9.2e} "
          f"{blkfail/N_FRAMES:8.3f} {rsfail/N_FRAMES:7.4f}")
