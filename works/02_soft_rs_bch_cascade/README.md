# Work 02 - Soft-Decision RS+BCH Cascade

## Objective

Explore a low-latency RS+BCH cascaded-code decoder over PAM4-AWGN.  The inner
BCH decoder is derived from Work 01; the novel system-level question is where
Lagrange interpolation can be shared across the outer RS and inner BCH stages.
The work compares a pure-RS baseline with cascade schemes A/B and LCC-BR/BM
outer decoding choices.

## Contents

| Location | Contents |
|---|---|
| `docs/advisor_brief/` | Research interpretation, work breakdown, and risks. |
| `docs/advisor_reply/` | Requirements clarification and design responses. |
| `docs/research_report/` | Problem definition and full feasibility/experimental report. |
| `code/python/cascade_src/` | Original PAM4, RS, Lagrange cache, cascade, and simulator code. |
| `code/python/experiments/` | n=63 smoke, n=255, all-config, KPI, and A/B-savings entry points. |
| `code/matlab/` | MATLAB reimplementation and tests. |
| `assets/data/` | Machine-readable benchmark and KPI results. |
| `assets/figures/` | FER/BER/KPI figures. |
| `assets/logs/` | Preserved original run logs. |

## Dependency boundary

This work imports the GF/BCH/LLOSD primitives from Work 01 explicitly through
repository-relative paths.  It does not mutate Work 01; all new simulation
results belong under this work's `assets/` directory.
