"""Import the explicitly shared BCH/OSD reproduction library.

Work 02 reuses the algorithms produced by Work 01.  The path is derived from
this file, so it remains valid regardless of the current working directory.
"""
import sys
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[5]
UPSTREAM_ROOT = REPOSITORY_ROOT / "works" / "01_bch_osd_reproduction" / "code" / "python"
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
