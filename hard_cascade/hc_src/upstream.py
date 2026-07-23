"""Re-export upstream GF + rs_bch_cascade's PureRSCodec / PAM4 modem."""
import sys
from pathlib import Path

# Upstream llosd_reproduction (has src/gf.py, src/bch.py, ...)
UPSTREAM = Path(__file__).parent.parent.parent
if str(UPSTREAM) not in sys.path:
    sys.path.insert(0, str(UPSTREAM))

# Also add rs_bch_cascade so we can import cascade_src.*
CASCADE = UPSTREAM / "rs_bch_cascade"
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
