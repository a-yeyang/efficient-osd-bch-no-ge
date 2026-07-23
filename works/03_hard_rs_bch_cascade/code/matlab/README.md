# Work 03 MATLAB reproduction

`HardCascade.m` is a Base MATLAB R2022b implementation of the Work 03 Python chain:

- primitive BCH(`2^m-1`, `2^m-1-2m`, `t=2`) systematic encoding;
- conventional (BM-equivalent locator + Chien) and direct LUT root-finding decoders;
- primitive systematic RS encoding and BM/Chien/Forney hard decoding;
- RS+BCH hard cascade and pure-RS baseline;
- LSB-first bit/symbol packing, Gray PAM4 hard decisions, Python-equivalent AWGN variance;
- operation counting, latency/KPI models, and Monte-Carlo simulation.

The code is self-contained within Work 03: it does not require Communications Toolbox, Python, or a neighboring work's MATLAB implementation.  Polynomial coefficients are low-degree first and systematic words are `[parity, message]`, matching the Python source exactly.

Run the real assertion suite from any current directory:

```matlab
addpath(genpath('.../works/03_hard_rs_bch_cascade/code/matlab'));
run_work03_tests
```

Experiment entries mirror the Python scripts:

```matlab
Work03Experiments.main(true);       % short profile; writes Work-03 assets
Work03Experiments.main(false);      % Python-sized Monte-Carlo schedule
Work03Experiments.main_v2(true);    % v2 latency experiment and figures
```

Pass `false` as the second argument to either entry point to suppress output files during automated tests.  Generated MATLAB data/figures use `matlab_profile_*` or `matlab_full_*` names under this work's `assets/` directory, so Python reference artifacts are never overwritten.

`reference_bch_t2/` remains a compact function-oriented BCH decoder reference. Its `smoke_test` now contains assertions for zero-, one-, and two-error correction plus a known Table-I uncorrectable branch.
