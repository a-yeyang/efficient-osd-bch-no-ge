"""Deterministic cross-work smoke tests for the preserved Python reference.

The paper-scale Monte-Carlo scripts intentionally remain separate: their
runtime is measured in hours and their random streams are not a stable
cross-language oracle.  This suite exercises the numerical paths that MATLAB
ports must match: finite-field/BCH operations, OSD-family decoders, PAM4/RS
cascades, BCH-t=2 direct decoding, and latency formulas.
"""
from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
WORK01 = REPOSITORY_ROOT / "works" / "01_bch_osd_reproduction" / "code" / "python"
WORK02 = REPOSITORY_ROOT / "works" / "02_soft_rs_bch_cascade" / "code" / "python"
WORK03 = REPOSITORY_ROOT / "works" / "03_hard_rs_bch_cascade" / "code" / "python"
for _path in (WORK03, WORK02, WORK01):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from src.bch import BCHCode, bpsk_modulate, llr_from_y  # noqa: E402
from src.decoders import hsd_fast, llosd_fast, sllosd_fast  # noqa: E402
from src.osd import osd_decode  # noqa: E402
from cascade_src.cascade import CascadeConfig, CascadedCodec, PureRSCodec  # noqa: E402
from cascade_src.pam4 import bits_to_pam4, pam4_to_bits_hard  # noqa: E402
from hc_src.bch_t2 import BCHt2Code  # noqa: E402
from hc_src.cascade import HardCascadeConfig, HardCascadedCodec, PureRSHardCodec  # noqa: E402
from hc_src.latency_model import LatencyModel  # noqa: E402


class PythonReferenceSmokeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.rng = np.random.default_rng(20260723)

    def test_work01_bch_and_osd_family(self) -> None:
        code = BCHCode(m=5, t=2)
        msg = self.rng.integers(0, 2, size=code.k, dtype=np.int8)
        codeword = code.encode(msg).astype(np.int8)
        self.assertTrue(np.all((codeword @ code.H.T) % 2 == 0))

        for positions in ([], [2], [1, 17]):
            received = codeword.copy()
            received[positions] ^= 1
            decoded, ok = code.bm_decode(received)
            self.assertTrue(ok)
            np.testing.assert_array_equal(decoded, codeword)

        llr = llr_from_y(bpsk_modulate(codeword), sigma=0.35)
        decoders = (
            lambda: osd_decode(code, llr, tau=1),
            lambda: llosd_fast(code, llr, tau=1),
            lambda: sllosd_fast(code, llr, theta_tuple=(2, 1)),
            lambda: hsd_fast(code, llr, tau=1, eta=1),
        )
        for decode in decoders:
            decoded, stats = decode()
            self.assertIsInstance(stats, dict)
            np.testing.assert_array_equal(decoded, codeword)
            self.assertTrue(np.all((decoded @ code.H.T) % 2 == 0))

    def test_work02_pam4_and_soft_cascade(self) -> None:
        bits = self.rng.integers(0, 2, size=32, dtype=np.int8)
        np.testing.assert_array_equal(pam4_to_bits_hard(bits_to_pam4(bits)), bits)

        cfg = CascadeConfig(m=5, k_rs=27, t_bch=1, llosd_tau=1, lcc_eta=1)
        message = self.rng.integers(0, 1 << cfg.m, size=cfg.k_rs, dtype=np.int64)
        cascade = CascadedCodec(cfg)
        encoded = cascade.encode(message)
        noiseless_llr = (1.0 - 2.0 * encoded) * 40.0
        for method in ("scheme_a", "scheme_b"):
            decoded, result = cascade.decode(noiseless_llr, method=method)
            self.assertTrue(result["ok"])
            np.testing.assert_array_equal(decoded, message)

        pure_rs = PureRSCodec(cfg)
        pure_encoded = pure_rs.encode(message)
        pure_llr = (1.0 - 2.0 * pure_encoded) * 40.0
        for method in ("hard", "soft"):
            decoded, result = pure_rs.decode(pure_llr, method=method)
            self.assertTrue(result["ok"])
            np.testing.assert_array_equal(decoded, message)

    def test_work03_direct_decoder_cascade_and_latency(self) -> None:
        bch = BCHt2Code(m=5)
        msg = self.rng.integers(0, 2, size=bch.k, dtype=np.int8)
        codeword = bch.encode(msg).astype(np.int8)
        for positions in ([], [3], [1, 13]):
            received = codeword.copy()
            received[positions] ^= 1
            conventional, conventional_ok = bch.decode_conventional(received)
            direct, direct_ok = bch.decode_direct(received)
            self.assertTrue(conventional_ok and direct_ok)
            np.testing.assert_array_equal(conventional, codeword)
            np.testing.assert_array_equal(direct, codeword)

        cfg = HardCascadeConfig(m=5, k_rs=27)
        message = self.rng.integers(0, 1 << cfg.m, size=cfg.k_rs, dtype=np.int64)
        for decoder in ("conv", "direct"):
            cascade = HardCascadedCodec(cfg, bch_decoder=decoder)
            decoded, ok, stats = cascade.decode(cascade.encode(message))
            self.assertTrue(ok)
            self.assertGreater(stats["n_bch_ok"], 0)
            np.testing.assert_array_equal(decoded, message)

        pure_rs = PureRSHardCodec(cfg)
        decoded, ok, _ = pure_rs.decode(pure_rs.encode(message))
        self.assertTrue(ok)
        np.testing.assert_array_equal(decoded, message)

        self.assertEqual(LatencyModel.bch_direct_cycles(2, m=8), 3)
        self.assertEqual(LatencyModel.bch_direct_cycles(2, m=7), 2)
        serial = LatencyModel.cascade_serial(3, cfg.t_rs)
        self.assertEqual(LatencyModel.cascade_lagrange_v1(3, cfg.t_rs), serial - 1)
        self.assertEqual(LatencyModel.cascade_lagrange_v2(3, cfg.t_rs), serial - 2)

    def test_relocated_experiment_sources_compile_and_paths_resolve(self) -> None:
        """Check every experiment source without triggering generation side effects.

        A few historical artifact scripts intentionally execute when invoked
        as scripts.  This test compiles all experiment sources and loads only
        their side-effect-free ``work_paths.py`` bootstrap modules, ensuring
        the reorganization is validated without overwriting reference assets.
        """
        experiment_dirs = (
            WORK01 / "experiments",
            WORK02 / "experiments",
            WORK03 / "experiments",
        )
        for index, experiment_dir in enumerate(experiment_dirs):
            for script in sorted(experiment_dir.glob("*.py")):
                if script.name == "work_paths.py":
                    continue
                compile(script.read_text(encoding="utf-8"), str(script), "exec")

            bootstrap = experiment_dir / "work_paths.py"
            spec = importlib.util.spec_from_file_location(
                f"work_paths_smoke_{index}", bootstrap
            )
            self.assertIsNotNone(spec and spec.loader)
            module = importlib.util.module_from_spec(spec)
            assert spec and spec.loader
            spec.loader.exec_module(module)
            expected_root = experiment_dir.parents[2]
            self.assertEqual(module.WORK_ROOT, expected_root)
            self.assertEqual(module.PYTHON_ROOT, expected_root / "code" / "python")
            self.assertEqual(module.ASSETS_ROOT, expected_root / "assets")
            self.assertTrue(module.ASSETS_ROOT.is_dir())


if __name__ == "__main__":
    unittest.main(verbosity=2)
