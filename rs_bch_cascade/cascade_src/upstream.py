"""Wrapper that re-exports the upstream (llosd_reproduction) codebase.

The upstream project at ~/workspace/llosd_reproduction uses `src/` as its
package name. We add its parent to sys.path so `import src.xxx` resolves to
the upstream code (there is no local `src/` in this subproject).
"""
import sys
from pathlib import Path

UPSTREAM_ROOT = Path(__file__).parent.parent.parent  # ~/workspace/llosd_reproduction
if str(UPSTREAM_ROOT) not in sys.path:
    sys.path.insert(0, str(UPSTREAM_ROOT))

from src.gf import GF, bch_generator_poly, bch_dimension, PRIM_POLY  # noqa: E402
from src.bch import BCHCode, bpsk_modulate, awgn, sigma_from_ebn0, llr_from_y  # noqa: E402
from src.osd import OpCounters  # noqa: E402
from src.decoders import llosd_fast  # noqa: E402
from src.llosd import build_rs_systematic_generator  # noqa: E402

__all__ = [
    "GF", "bch_generator_poly", "bch_dimension", "PRIM_POLY",
    "BCHCode", "bpsk_modulate", "awgn", "sigma_from_ebn0", "llr_from_y",
    "OpCounters", "llosd_fast", "build_rs_systematic_generator",
]
