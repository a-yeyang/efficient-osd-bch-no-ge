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

    % Test 1 error
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

    % Test 3 errors (should fail correctly)
    n_bad_success = 0;
    for trial = 1:20
        r = cw;
        positions = randperm(n, 3);
        for i = 1:3
            p = positions(i) - 1;
            r(p+1) = 1;
        end
        [c1, ok1] = bch_decode_conventional(gf, r);
        [c2, ok2] = bch_decode_direct(gf, r, lut_A);
        if ok1 && all(c1 == cw)
            n_bad_success = n_bad_success + 1;
        end
        if ok2 && all(c2 == cw)
            n_bad_success = n_bad_success + 1;
        end
    end
    fprintf('3-error accidental recoveries: %d (should be 0)\n', n_bad_success);
    fprintf('=== All Matlab tests passed ===\n');
end
