# Work 01 - BCH OSD Paper Reproduction

## Objective

Reproduce the IEEE TIT 2025 paper *Efficient Ordered Statistics Decoding of
BCH Codes Without Gaussian Elimination*.  The key idea is to avoid the
Gaussian-elimination bottleneck in ordered-statistics decoding by using a
Lagrange-interpolated systematic Reed-Solomon generator representation.  The
work implements LLOSD, LLOSD-B, SLLOSD, and HSD alongside OSD/YSVL/CJ/PLCC
baselines.

## Contents

| Location | Contents |
|---|---|
| `docs/reference_paper/` | The cited IEEE paper. |
| `docs/reproduction_report.md` | Chinese analysis, result comparison, known deviations, and methodology. |
| `docs/latency_analysis/` | Latency analysis source and rendered report. |
| `docs/*.docx` | Generated reproduction report. |
| `code/python/src/` | Original Python GF, BCH, OSD-family, and Monte-Carlo implementation. |
| `code/python/experiments/` | One reproducibility entry point per figure/table/report artifact. |
| `code/matlab/` | MATLAB reimplementation and deterministic/smoke tests. |
| `assets/data/` | Figure data in JSON. |
| `assets/figures/` | Reproduced PNG/PDF figures. |
| `assets/tables/` | Table I-IV JSON and Markdown results. |

## Outputs and traceability

`fig02_nbch.py` produces Fig. 2; `fig03_04_fer.py` produces Figs. 3-4;
`fig05_comparison.py` produces Fig. 5; `fig07_09_longcode.py` and
`fig07_09_hsd_rerun.py` produce Figs. 7 and 9; `fig08_rate_sweep.py` produces
Fig. 8; `fig10_11_ops.py` produces Figs. 10-11; and `tables.py` produces
Tables I-IV.  `latency_analysis.py` reads the table JSON to create the latency
figures.

Python output paths are resolved from the experiment file and write only to
this work's `assets/` directory.  The MATLAB test runner is under
`code/matlab/tests/`.
