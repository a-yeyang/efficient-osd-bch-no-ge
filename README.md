# BCH/RS-BCH Reproduction and Research Workspace

This repository is organized as a sequence of three linked research works.  Each
work owns its source code, papers/reports, experimental data, and figures.  The
dependency direction is intentional: later work may reuse an earlier work, but
does not modify it.

```text
Work 01: Efficient OSD of BCH codes without Gaussian elimination
    └── supplies GF/BCH/OSD/LLOSD building blocks
          ↓
Work 02: Soft-decision RS+BCH cascade with Lagrange sharing
    └── reuses Work 01 and adds PAM4, RS, and soft cascade simulation
          ↓
Work 03: Hard-decision RS+BCH cascade low-latency decoder
    └── reuses Work 01/02 concepts and adds BCH-t=2 direct decoding
        and hardware-oriented latency modelling
```

## Directory map

| Work | Scope | Primary materials |
|---|---|---|
| [`works/01_bch_osd_reproduction`](works/01_bch_osd_reproduction) | IEEE TIT 2025 BCH-OSD paper reproduction | `code/python`, `code/matlab`, `assets/{data,figures,tables}`, `docs/` |
| [`works/02_soft_rs_bch_cascade`](works/02_soft_rs_bch_cascade) | Soft-decision RS+BCH cascade research | `code/python`, `code/matlab`, `assets/{data,figures,logs}`, `docs/` |
| [`works/03_hard_rs_bch_cascade`](works/03_hard_rs_bch_cascade) | Hard-decision RS+BCH cascade and latency research | `code/python`, `code/matlab`, `assets/{data,figures,logs}`, `docs/` |

Within a work:

- `code/python/` preserves the original Python implementation and its experiments.
- `code/matlab/` contains the MATLAB reimplementation and MATLAB tests.
- `assets/` contains reproducibility artifacts, not source code.
- `docs/` contains the cited paper, analysis reports, advisor correspondence, and feasibility material.

## Reproducibility entry points

Python and MATLAB paths are derived from the script locations, so commands can
be launched from any current directory.  Use the bundled or a local Python
environment with the dependencies listed in `requirements.txt`.

MATLAB R2022b or later is supported.  Each work has a test runner under
`code/matlab/tests/`; the repository-level runner is `run_all_matlab_tests.m`.
It runs deterministic unit tests and small Monte-Carlo smoke tests.  The full
paper-scale Monte-Carlo configurations remain available in the experiment
entry points and deliberately are not run by the smoke suite.

From the repository root, run the delivered acceptance checks with:

```powershell
.\.venv\Scripts\python.exe -m unittest discover -s tests -v
matlab -batch "run('run_all_matlab_tests.m')"
```

## Research provenance

Work 01 reproduces the paper *Efficient Ordered Statistics Decoding of BCH
Codes Without Gaussian Elimination*.  Work 02 investigates Lagrange sharing
in a soft-decision PAM4 RS+BCH cascade.  Work 03 follows with a hard-decision,
direct-root-finding BCH-t=2 design and a cycle-level latency model.  See each
work's `docs/` directory for the underlying paper/report and the assumptions
behind reported figures and KPIs.

## License

The implementation is released under the repository's MIT license.  The
included IEEE paper and other authored reports retain their respective rights.
