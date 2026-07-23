% Smoke test for Matlab BCH decoders. Run in Matlab / Octave:
%   >> smoke_test
%
% Expected output: all 0-error / 1-error / 2-error tests pass.

function smoke_test()
    fprintf('=== Matlab BCH t=2 Decoder Smoke Test ===\n');
    gf = GF_init(8);
    lut_A = build_lut_A(gf);
    fprintf('GF(2^%d), n=%d, LUT_A built (valid entries: %d/%d)\n', ...
            gf.m, gf.n, sum(lut_A.valid), length(lut_A.valid));

    n = gf.n;
    % Build a valid BCH codeword by hand: c = 0 (all-zeros) satisfies S_i = 0.
    rng(42);
    cw = zeros(1, n);

    % Test 0 errors
    r = cw;
    [c1, ok1] = bch_decode_conventional(gf, r);
    [c2, ok2] = bch_decode_direct(gf, r, lut_A);
    assert(ok1 && ok2 && all(c1 == cw) && all(c2 == cw), 'Zero-error failed');
    fprintf('0-error: OK\n');

    % Test 1 error.  This is an assertion-based smoke test: a printed count
    % is not treated as success unless every deterministic pattern passed.
    n_pass = 0;
    for trial = 1:20
        r = cw;
        p = randi(n) - 1;
        r(p+1) = 1;
        [c1, ok1] = bch_decode_conventional(gf, r);
        [c2, ok2] = bch_decode_direct(gf, r, lut_A);
        if ok1 && ok2 && all(c1 == cw) && all(c2 == cw)
            n_pass = n_pass + 1;
        end
    end
    fprintf('1-error: %d/20 pass\n', n_pass);
    assert(n_pass == 20, 'One-error correction failed for at least one decoder/pattern.');

    % Test 2 errors
    n_pass = 0;
    for trial = 1:30
        r = cw;
        positions = randperm(n, 2);
        for i = 1:2
            p = positions(i) - 1;
            r(p+1) = 1;
        end
        [c1, ok1] = bch_decode_conventional(gf, r);
        [c2, ok2] = bch_decode_direct(gf, r, lut_A);
        if ok1 && ok2 && all(c1 == cw) && all(c2 == cw)
            n_pass = n_pass + 1;
        end
    end
    fprintf('2-error: %d/30 pass\n', n_pass);
    assert(n_pass == 30, 'Two-error correction failed for at least one decoder/pattern.');

    % A generic three-error pattern can be miscorrected by any bounded-
    % distance decoder.  Instead, construct the specific Table-I failure
    % class with error locators {1, alpha, 1+alpha}; then S1=0 and S3~=0,
    % so both decoders must explicitly report an uncorrectable word.
    alpha = gf.EXP(2);   % alpha^1 (MATLAB's one-based EXP indexing)
    third_locator = bitxor(1, alpha);
    positions = [0, 1, gf.LOG(third_locator + 1)];
    r = cw;
    r(positions + 1) = 1;
    [~, ok1] = bch_decode_conventional(gf, r);
    [~, ok2] = bch_decode_direct(gf, r, lut_A);
    assert(~ok1 && ~ok2, 'Known S1=0/S3~=0 three-error branch was accepted.');
    fprintf('Known uncorrectable 3-error pattern: OK\n');
    fprintf('=== All Matlab tests passed ===\n');
end
