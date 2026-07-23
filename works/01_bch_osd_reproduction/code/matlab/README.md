# Work 01 MATLAB reproduction

The canonical Base MATLAB R2022b port is the `+work01` package:

| Python reference | Canonical MATLAB counterpart |
|---|---|
| `src/gf.py` | `work01.Core.gf_*` finite-field and polynomial methods |
| `src/bch.py` | `work01.Core.bch_*`, BPSK/AWGN helpers, and BM decoder |
| `src/osd.py` | `work01.Core.gaussian_elim_binary` and `osd_decode` |
| `src/llosd*.py`, `src/decoders.py` | `build_rs_systematic_generator`, `llosd_decode`, and `llosd_fast` |
| `src/sllosd*.py` | `sllosd_fast` and theoretical candidate counts |
| `src/hsd.py`, `src/baselines.py`, `src/ml.py`, `src/simulate.py` | `hsd_fast`, baseline wrappers, ML search, and `run_mc` |
| Python figure/table scripts | matching functions in `experiments/`, backed by `work01.Experiment` |
| Python DOCX helpers | `experiments/build_docx_report.m` and `embed_pdf_in_docx.m` (Windows + Word) |

Run the deterministic acceptance suite from the repository root:

```matlab
run('run_all_matlab_tests.m')
```

Or test this work alone:

```matlab
addpath(genpath('works/01_bch_osd_reproduction/code/matlab'));
run_work01_tests
```

The experiment functions retain the Python-sized Monte-Carlo schedules; the
acceptance suite deliberately uses fixed vectors and tiny simulations.  The
generated Word artifacts are written to `../../assets/generated/` and never
overwrite the original report under `../../docs/`.

`BCHOSD.m` and `Work01Experiments.m` are retained alternate prototypes from
the working tree.  They are not the canonical port and are not used by
`run_work01_tests` or the repository-level runner.
