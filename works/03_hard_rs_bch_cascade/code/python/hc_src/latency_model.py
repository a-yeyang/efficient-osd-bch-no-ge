"""Latency model with two metrics: F_{2^m} ops and clock cycles.

Clock cycle model follows Lagendijk 2026 Table VI for n=256:
    BCH Conv (BM+Chien) t=2:  8 cycles
    BCH Direct t=2:            3 cycles
    RS-BM t=8:                 depends on t; roughly 2t + 2 for BM iterations
                               + 1 for syndrome + 1 for Chien + 1 for Forney

For fairness we model:
    RS-BM(n, t) = 1 (syndrome) + 2t (BM iterations) + 1 (Chien) + 1 (Forney)
                = 2t + 3 cycles

    BCH-Conv(t=2) = 8  (paper value for n=256)
    BCH-Direct(t=2) at n=255 (GF(2^8)) = 3 (paper value)
    BCH-Direct(t=2) at n=127 (GF(2^7)) = 2 (§6.3 optimization: smaller LUT, deeper
                                            combinational logic can fit into 2 cycles)

We assume all operations within a cycle are fully parallel.

Lagrange sharing modes:
    'none':      no sharing. Cascade = T_BCH + T_RS_BM
    'v1':        share alpha table + syndrome unit. Save 1 cycle.
                 Cascade = T_BCH + T_RS_BM - 1
    'v2':        v1 + share GF multiplier array (BCH LUT + RS Chien can time-
                 multiplex the same hardware). Save 1 additional cycle.
                 Cascade = T_BCH + T_RS_BM - 2
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
    def bch_conv_cycles(t: int, m: int = 8) -> int:
        """BCH Conventional BM + Chien for binary BCH."""
        if t == 2:
            return 8
        elif t == 3:
            return 8
        elif t == 4:
            return 10
        return 2 * t + 2 + 2  # extrapolate

    @staticmethod
    def bch_direct_cycles(t: int, m: int = 8) -> int:
        """BCH Direct root finding.

        For n=256 (GF(2^8)) the paper reports 3 cycles.
        For n=127 (GF(2^7)) the smaller LUT (128 vs 256 entries) has shallower
        combinational depth; precompute + LUT lookup can be merged into 1 cycle,
        yielding 2 cycles total (§6.3 Multi-cycle refinement).
        """
        if t == 2:
            if m <= 7:
                return 2  # optimized for GF(2^7) or smaller
            else:
                return 3  # paper value for GF(2^8)
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
    def cascade_lagrange_v1(cls, bch_cycles: int, rs_t: int) -> int:
        """Cascade with hardware Lagrange sharing v1.

        Share the α power LUT + syndrome computation unit.
        The syndrome cycle of RS-BM can overlap with the last cycle of BCH.
        Save 1 cycle.
        """
        return bch_cycles + cls.rs_bm_cycles(rs_t) - 1

    @classmethod
    def cascade_lagrange_v2(cls, bch_cycles: int, rs_t: int) -> int:
        """Cascade with aggressive Lagrange sharing v2.

        On top of v1, additionally share the GF multiplier array:
          - BCH Direct's LUT-lookup + correction stage and
            RS-BM's Chien search stage
          both need parallel GF multipliers. They can time-multiplex the same
          hardware in adjacent cycles.
        Save 1 additional cycle → total 2 cycles saved vs 'none'.
        """
        return bch_cycles + cls.rs_bm_cycles(rs_t) - 2

    @classmethod
    def cascade_lagrange_shared(cls, bch_cycles: int, rs_t: int) -> int:
        """Backward-compat alias for v1."""
        return cls.cascade_lagrange_v1(bch_cycles, rs_t)

    # -----------------------------------------------------------------
    # KPI check
    # -----------------------------------------------------------------
    @classmethod
    def kpi_check(cls, cascade_cycles: int, baseline_cycles: int) -> tuple:
        """Return (ratio_pct, passes_10pct).
        ratio_pct = (cascade - baseline) / baseline * 100."""
        ratio = (cascade_cycles - baseline_cycles) / baseline_cycles * 100
        return ratio, ratio <= 10.0
