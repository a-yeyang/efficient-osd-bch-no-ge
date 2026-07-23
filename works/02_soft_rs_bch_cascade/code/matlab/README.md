# Work 02 MATLAB reproduction

This directory is a Base MATLAB R2022b implementation of the Work 02 Python
reference. Start MATLAB with this folder on its path, then call
`startup_work02`.

Core APIs are in the `+cascade` package: PAM4/AWGN/LLR, GF/BCH/LLOSD, RS
BM/LCC-BR, Scheme A/B cascade decoding, Lagrange caching, and Monte-Carlo
simulation. The five experiment entry points are in `experiments/`:

- `run_smoke_n63('full'|'quick')`
- `run_n255_scheme_a('full'|'quick')`
- `run_all_configs('full'|'quick')`
- `run_scheme_ab_savings('full'|'quick')`
- `run_kpi_analysis('full'|'quick')`

All entry points derive output paths from their own location and write only to
Work 02 `assets/data` and `assets/figures`. Run `run_work02_tests` for the
fast assertion suite.

Scheme A's operation counter includes every inner BCH/LLOSD invocation as
well as outer RS work. Scheme B uses one `LagrangeCache` through its inner
Lagrange generator and outer RS locator evaluations; its reported saving is
derived from observed denominator-cache reuse, with cache statistics returned
in the decoder result. PAM4 uses the same reference convention as Python:
`sigma^2 = 5 * rate / (4 * 10^(EbN0_dB/10))`.

`SoftCascade.m` and `Work02Experiments.m` are retained alternate prototypes.
The canonical, tested API is the `+cascade` package and the `experiments/`
entry points above; the repository-level runner does not use the prototypes.
