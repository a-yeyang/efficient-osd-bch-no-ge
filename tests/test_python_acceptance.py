"""Deterministic post-reorganization acceptance tests for all three works."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "works" / "01_bch_osd_reproduction" / "code" / "python"))
sys.path.insert(0, str(ROOT / "works" / "02_soft_rs_bch_cascade" / "code" / "python"))
sys.path.insert(0, str(ROOT / "works" / "03_hard_rs_bch_cascade" / "code" / "python"))

from src.bch import BCHCode
from src.decoders import hsd_fast, llosd_fast, sllosd_fast
from src.osd import osd_decode
from cascade_src.cascade import CascadeConfig, CascadedCodec
from hc_src.bch_t2 import BCHt2Code
from hc_src.cascade import HardCascadeConfig, HardCascadedCodec


class AcceptanceTests(unittest.TestCase):
    def test_work01_decoders_and_bounded_distance(self) -> None:
        code = BCHCode(5, 2)
        msg = np.arange(code.k, dtype=np.int8) % 2
        word = code.encode(msg)
        for positions in [[], [0], [0, 7], [4, 21]]:
            received = word.copy()
            received[positions] ^= 1
            decoded, ok = code.bm_decode(received)
            self.assertTrue(ok)
            np.testing.assert_array_equal(decoded, word)
        llr = (1.0 - 2.0 * word) * 12.0
        for fn, args in [
            (osd_decode, (1,)),
            (llosd_fast, (1,)),
            (sllosd_fast, ((2, 1),)),
            (hsd_fast, (1, 3)),
        ]:
            decoded, _ = fn(code, llr, *args)
            np.testing.assert_array_equal(decoded, word, err_msg=fn.__name__)

    def test_work02_soft_cascade_noiseless(self) -> None:
        cfg = CascadeConfig(m=5, k_rs=27, t_bch=1, llosd_tau=1, lcc_eta=2)
        codec = CascadedCodec(cfg)
        msg = np.arange(cfg.k_rs, dtype=np.int64) % (1 << cfg.m)
        bits = codec.encode(msg)
        decoded, result = codec.decode((1.0 - 2.0 * bits) * 30.0, "scheme_a")
        self.assertTrue(result["ok"])
        np.testing.assert_array_equal(decoded, msg)

    def test_work03_direct_and_hard_cascade(self) -> None:
        bch = BCHt2Code(5)
        msg = np.arange(bch.k, dtype=np.int8) % 2
        word = bch.encode(msg)
        for decoder in (bch.decode_conventional, bch.decode_direct):
            decoded, ok = decoder(word.copy())
            self.assertTrue(ok)
            np.testing.assert_array_equal(decoded, word)
        errored = word.copy()
        errored[[1, 6]] ^= 1
        decoded, ok = bch.decode_direct(errored)
        self.assertTrue(ok)
        np.testing.assert_array_equal(decoded, word)
        cfg = HardCascadeConfig(m=5, k_rs=27)
        codec = HardCascadedCodec(cfg, "direct")
        symbols = np.arange(cfg.k_rs, dtype=np.int64) % (1 << cfg.m)
        decoded, ok, _ = codec.decode(codec.encode(symbols))
        self.assertTrue(ok)
        np.testing.assert_array_equal(decoded, symbols)


if __name__ == "__main__":
    unittest.main(verbosity=2)
