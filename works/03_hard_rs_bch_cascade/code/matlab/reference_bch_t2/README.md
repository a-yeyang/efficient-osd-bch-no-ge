# BCH t=2 function-oriented reference

This folder preserves a compact, standalone function-oriented implementation of the two BCH t=2 hard-decision decoders described in Work 03:

- `bch_decode_conventional.m`: the conventional BM-equivalent locator plus Chien search;
- `bch_decode_direct.m`: direct root finding through the `A(Y)=Y^2+Y+k` lookup table;
- `GF_init.m`, `gf_mul.m`, and `gf_div.m`: GF(2^m) EXP/LOG arithmetic;
- `bch_syndromes.m` and `build_lut_A.m`: shared odd-syndrome/LUT machinery.

Run the focused assertion smoke test with:

```matlab
addpath('.../works/03_hard_rs_bch_cascade/code/matlab/reference_bch_t2');
smoke_test
```

It asserts zero-, one-, and two-error correction for both decoders, and tests a deliberately constructed `S1=0, S3~=0` uncorrectable three-error branch.  The complete self-contained RS/BCH cascade, PAM4 channel, simulation driver, latency models, and full tests live one level above in `HardCascade.m`, `Work03Experiments.m`, and `tests/run_work03_tests.m`.
