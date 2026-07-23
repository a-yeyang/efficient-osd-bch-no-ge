# Work 03 - Hard-Decision RS+BCH Cascade

## Objective

Develop a hardware-oriented hard-decision RS+BCH cascade.  It replaces the
soft inner decoder with a BCH-t=2 direct root-finding decoder and evaluates
cycle-level latency.  The v2 report further investigates aggressive Lagrange
sharing and a refined two-cycle direct decoder for n=127.

## Contents

| Location | Contents |
|---|---|
| `docs/advisor_reply/` | Lagendijk direct-decoding interpretation and requirement analysis. |
| `docs/research_report/` | Original and v2 feasibility/implementation reports. |
| `code/python/hc_src/` | Original BCH-t=2, hard cascade, simulator, and latency model. |
| `code/python/experiments/` | Original FER/KPI and v2 latency experiment drivers. |
| `code/matlab/reference_bch_t2/` | Preserved initial MATLAB BCH-t=2 reference port. |
| `code/matlab/` | Complete MATLAB reimplementation and tests. |
| `assets/data/` | Benchmark and v2 latency data. |
| `assets/figures/` | FER and latency figures. |
| `assets/logs/` | Preserved original run log. |

## Dependency boundary

This work reuses Work 01's GF/BCH foundations and Work 02's RS/PAM4 concepts
through repository-relative imports.  Its direct decoder and latency model are
kept local because they are the research contribution of this phase.
