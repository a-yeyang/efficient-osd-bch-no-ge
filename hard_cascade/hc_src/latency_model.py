"""Latency model with two metrics: F_{2^m} ops and clock cycles.

Clock cycle model follows Lagendijk 2026 Table VI for n=256:
    BCH Conv (BM+Chien) t=2:  8 cycles
    BCH Direct t=2:            3 cycles
    RS-BM t=8:                 depends on t; roughly 2t + 2 for BM iterations
                               + 1 for syndrome + 1 for Chien + 1 for Forney

For fairness we model:
    RS-BM(n, t) = 1 (syndrome) + 2t (BM iterations) + 1 (Chien) + 1 (Forney)
                = 2t + 3 cycles

    BCH-Conv(t=2) = 8  (paper value)
    BCH-Direct(t=2) = 3 (paper value)

We assume all operations within a cycle are fully parallel (P → ∞ semantics
for that cycle's operations, so cycle count is decoupled from ops count).

For hardware Lagrange sharing (Q5): when inner + outer share the α power
LUT and one syndrome evaluation unit, the syndrome computation cycle
overlaps (both compute their syndromes in parallel in the same cycle).
This saves ~1 cycle on the cascade timeline.
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass
class LatencyModel:
    """Clock cycle model for RS-BM and BCH decoders (Lagendijk-style)."""

    @staticmethod
    def rs_bm_cycles(t: int) -> int:
        """RS BM + Chien + Forney: 1 (syndrome) + 2t (BM) + 1 (Chien) + 1 (Forney) = 2t+3."""
        return 2 * t + 3

    @staticmethod
    def bch_conv_cycles(t: int) -> int:
        """BCH Conventional BM + Chien for binary BCH."""
        if t == 2:
            return 8
        elif t == 3:
            return 8
        elif t == 4:
            return 10
        return 2 * t + 2 + 2  # extrapolate

    @staticmethod
    def bch_direct_cycles(t: int) -> int:
        """BCH Direct root finding."""
        if t == 2:
            return 3
        elif t == 3:
            return 4
        elif t == 4:
            return 8
        return None

    # -----------------------------------------------------------------
    # Cascade timing (serial pipeline: BCH inner → RS outer)
    # -----------------------------------------------------------------
    @classmethod
    def cascade_serial(cls, bch_cycles: int, rs_t: int) -> int:
        """Cascade with serial pipeline: BCH inner then RS outer."""
        return bch_cycles + cls.rs_bm_cycles(rs_t)

    @classmethod
    def cascade_lagrange_shared(cls, bch_cycles: int, rs_t: int) -> int:
        """Cascade with hardware Lagrange sharing.

        Assumption: BCH and RS share the α power LUT + syndrome unit,
        so the syndrome computation cycle can be overlapped. Save 1 cycle.
        """
        return bch_cycles + cls.rs_bm_cycles(rs_t) - 1

    # -----------------------------------------------------------------
    # KPI check
    # -----------------------------------------------------------------
    @classmethod
    def kpi_check(cls, cascade_cycles: int, baseline_cycles: int) -> tuple:
        """Return (ratio, passes_10pct).
        ratio = (cascade - baseline) / baseline * 100 (%)."""
        ratio = (cascade_cycles - baseline_cycles) / baseline_cycles * 100
        return ratio, ratio <= 10.0
