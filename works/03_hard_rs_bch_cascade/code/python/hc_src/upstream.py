"""Import shared Work 01 algebra and Work 02 RS/PAM4 components explicitly."""
import sys
from pathlib import Path

# Reuse the prior works through stable repository-relative paths.
REPOSITORY_ROOT = Path(__file__).resolve().parents[5]
UPSTREAM = REPOSITORY_ROOT / "works" / "01_bch_osd_reproduction" / "code" / "python"
if str(UPSTREAM) not in sys.path:
    sys.path.insert(0, str(UPSTREAM))

# Work 02 owns the RS and PAM4 implementation.
CASCADE = REPOSITORY_ROOT / "works" / "02_soft_rs_bch_cascade" / "code" / "python"
if str(CASCADE) not in sys.path:
    sys.path.insert(0, str(CASCADE))

from src.gf import GF, bch_generator_poly, bch_dimension, PRIM_POLY  # noqa: E402
from src.bch import bpsk_modulate, awgn, sigma_from_ebn0, llr_from_y  # noqa: E402
from src.osd import OpCounters  # noqa: E402

from cascade_src.pam4 import (  # noqa: E402
    bits_to_pam4, pam4_to_bits_hard, sigma_from_ebn0_pam4,
    awgn_channel, pam4_bit_llr, E_S_AVG, E_B_AVG,
)
from cascade_src.rs_code import RSCode  # noqa: E402
from cascade_src.cascade import run_channel  # noqa: E402

__all__ = [
    "GF", "bch_generator_poly", "bch_dimension", "PRIM_POLY",
    "bpsk_modulate", "awgn", "sigma_from_ebn0", "llr_from_y",
    "OpCounters",
    "bits_to_pam4", "pam4_to_bits_hard", "sigma_from_ebn0_pam4",
    "awgn_channel", "pam4_bit_llr", "E_S_AVG", "E_B_AVG",
    "RSCode", "run_channel",
]
