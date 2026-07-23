"""Stable paths for Work 01 experiment entry points."""
from pathlib import Path
import os
import sys

WORK_ROOT = Path(__file__).resolve().parents[3]
PYTHON_ROOT = WORK_ROOT / "code" / "python"
ASSETS_ROOT = WORK_ROOT / "assets"
DATA_DIR = ASSETS_ROOT / "data"
FIGURES_DIR = ASSETS_ROOT / "figures"
TABLES_DIR = ASSETS_ROOT / "tables"


def setup() -> None:
    """Make legacy relative experiment outputs resolve under ``assets/``."""
    if str(PYTHON_ROOT) not in sys.path:
        sys.path.insert(0, str(PYTHON_ROOT))
    os.chdir(ASSETS_ROOT)
