function run_work03_tests()
%RUN_WORK03_TESTS Deterministic, assertion-based MATLAB validation for Work 03.
%
% This runner intentionally uses only paths derived from its own location.  It
% validates nonzero BCH words, every <=2-error pattern for BCH(31,21), direct
% and conventional decoders, RS/pure-cascade paths, PAM4, both latency models,
% the bounded main profile, and v2 artifact output.

    thisFile = mfilename('fullpath');
    matlabRoot = fileparts(fileparts(thisFile));
    addpath(genpath(matlabRoot));
    rng(20260723, 'twister');

    fprintf('Work 03 MATLAB tests: BCH(t=2), RS cascade, PAM4, latency, experiments\n');

    %% BCH(31,21), nonzero systematic word and complete <=2-error coverage.
    bch = HardCascade.bch_t2_create(5, true);
    assert(bch.n == 31 && bch.k == 21 && bch.t == 2, 'Unexpected BCH(31,21) parameters.');
    assert(sum(bch.lut.valid) == 16, 'GF(2^5) direct-root LUT should have 16 valid entries.');
    message = mod(0:bch.k-1, 2);
    word = HardCascade.bch_t2_encode(bch, message);
    assert(any(word), 'The BCH test word must be nonzero.');
    assert(isequal(word, mod(message * bch.G, 2)), 'BCH polynomial and G-matrix encoders disagree.');
    [s1, s3] = HardCascade.bch_syndromes(bch, word);
    assert(s1 == 0 && s3 == 0, 'Encoded nonzero BCH word has nonzero syndrome.');

    decoderFns = {@HardCascade.bch_decode_conventional, @HardCascade.bch_decode_direct};
    decoderNames = {'conventional', 'direct'};
    for decoderIndex = 1:numel(decoderFns)
        [decoded, ok, counter] = decoderFns{decoderIndex}(bch, word);
        assert(ok && isequal(decoded, word) && counter.f2m >= 0, ...
            '%s decoder failed the zero-error word.', decoderNames{decoderIndex});
        for first = 1:bch.n
            received = word;
            received(first) = 1 - received(first);
            [decoded, ok] = decoderFns{decoderIndex}(bch, received);
            assert(ok && isequal(decoded, word), '%s decoder failed a one-error word.', decoderNames{decoderIndex});
        end
        for first = 1:bch.n-1
            for second = first+1:bch.n
                received = word;
                received([first second]) = 1 - received([first second]);
                [decoded, ok] = decoderFns{decoderIndex}(bch, received);
                assert(ok && isequal(decoded, word), ...
                    '%s decoder failed pair (%d,%d).', decoderNames{decoderIndex}, first, second);
            end
        end
    end

    % Exercise the actual n=127/n=255 BCH dimensions with nonzero words.
    for m = [7 8]
        longBch = HardCascade.bch_t2_create(m, true);
        longMessage = mod(0:longBch.k-1, 2);
        longWord = HardCascade.bch_t2_encode(longBch, longMessage);
        assert(any(longWord), 'Long BCH test word unexpectedly zero.');
        for positions = {[1], [2, longBch.n]}
            received = longWord;
            selected = positions{1};
            received(selected) = 1 - received(selected);
            [directWord, directOk] = HardCascade.bch_decode_direct(longBch, received);
            [conventionalWord, conventionalOk] = HardCascade.bch_decode_conventional(longBch, received);
            assert(directOk && conventionalOk && isequal(directWord, longWord) && ...
                isequal(conventionalWord, longWord), 'Long BCH bounded-distance decode failed.');
        end
    end

    % Deliberately construct the Table-I S1=0, S3~=0 uncorrectable branch.
    alpha = bch.gf.EXP(2);
    locator3 = bitxor(1, alpha);
    positions = [0, 1, bch.gf.LOG(locator3 + 1)] + 1;
    received = word;
    received(positions) = 1 - received(positions);
    [~, conventionalOk] = HardCascade.bch_decode_conventional(bch, received);
    [~, directOk] = HardCascade.bch_decode_direct(bch, received);
    assert(~conventionalOk && ~directOk, 'Known uncorrectable BCH branch was accepted.');

    %% RS BM and pure-RS hard baseline.
    rs = HardCascade.rs_create(5, 27);
    rsMessage = mod(0:rs.k-1, 2^rs.m);
    rsWord = HardCascade.rs_encode_systematic(rs, rsMessage);
    [decodedRs, rsOk] = HardCascade.rs_bm_decode(rs, rsWord);
    assert(rsOk && isequal(decodedRs, rsWord), 'Noiseless RS BM decode failed.');
    erroredRs = rsWord;
    erroredRs([1 8]) = bitxor(erroredRs([1 8]), [5 11]);
    [decodedRs, rsOk] = HardCascade.rs_bm_decode(rs, erroredRs);
    assert(rsOk && isequal(decodedRs, rsWord), 'RS BM failed to correct two symbol errors.');

    pure = HardCascade.pure_rs_create(5, 27);
    pureBits = HardCascade.pure_rs_encode(pure, rsMessage);
    % One bit in each of two distinct RS symbols creates two symbol errors.
    pureDamaged = pureBits;
    pureDamaged([1, 5*9 + 1]) = 1 - pureDamaged([1, 5*9 + 1]);
    [pureMessage, pureOk] = HardCascade.pure_rs_decode(pure, pureDamaged);
    assert(pureOk && isequal(pureMessage, rsMessage), 'Pure-RS hard baseline failed bounded-distance decoding.');

    %% Full direct/conventional cascade, including two errors in every BCH block.
    for decoder = {'direct', 'conv'}
        cascade = HardCascade.cascade_create(5, 27, decoder{1});
        coded = HardCascade.cascade_encode(cascade, rsMessage);
        [decodedMessage, cascadeOk, stats] = HardCascade.cascade_decode(cascade, coded);
        assert(cascadeOk && isequal(decodedMessage, rsMessage) && ...
            stats.n_bch_ok == cascade.n_bch_blocks, 'Noiseless %s cascade failed.', decoder{1});
        damaged = coded;
        for block = 1:cascade.n_bch_blocks
            indices = (block - 1) * cascade.bch.n + [2 7];
            damaged(indices) = 1 - damaged(indices);
        end
        [decodedMessage, cascadeOk, stats] = HardCascade.cascade_decode(cascade, damaged);
        assert(cascadeOk && isequal(decodedMessage, rsMessage) && ...
            stats.n_bch_ok == cascade.n_bch_blocks, ...
            '%s cascade failed two correctable errors per inner block.', decoder{1});
    end

    %% PAM4 map/channel and Python-equivalent sigma law.
    pamBits = [0 0 0 1 1 1 1 0];
    assert(isequal(HardCascade.bits_to_pam4(pamBits), [-3 -1 1 3]), 'PAM4 mapper mismatch.');
    assert(isequal(HardCascade.pam4_to_bits_hard([-3.2 -0.8 0.8 3.2]), pamBits), ...
        'PAM4 nearest-level hard demapper mismatch.');
    assert(isequal(HardCascade.run_pam4_channel_hard(pamBits, 100, 0.5, 42), pamBits), ...
        'High-SNR hard PAM4 channel changed a deterministic word.');
    assert(abs(HardCascade.sigma_from_ebn0_pam4(0, 0.5) - sqrt(5/8)) < 1e-12, ...
        'PAM4 sigma does not match sigma^2=5*R/(4*EbN0).');

    %% Exact latency-model values from main.py and main_v2.py.
    cfg255 = HardCascade.hard_cascade_config(8, 239);
    cfg127 = HardCascade.hard_cascade_config(7, 113);
    v1_255 = HardCascade.latency_summary_v1(cfg255);
    v1_127 = HardCascade.latency_summary_v1(cfg127);
    v2_255 = HardCascade.latency_summary_v2(cfg255);
    v2_127 = HardCascade.latency_summary_v2(cfg127);
    assert(isequal(v1_255, struct('pure_rs', 19, 'conv_no_share', 27, ...
        'direct_no_share', 22, 'conv_lagrange', 26, 'direct_lagrange', 21, 'kpi_target', 20)), ...
        'v1 n255 latency mismatch.');
    assert(isequal(v1_127, struct('pure_rs', 17, 'conv_no_share', 25, ...
        'direct_no_share', 20, 'conv_lagrange', 24, 'direct_lagrange', 19, 'kpi_target', 18)), ...
        'v1 n127 latency mismatch.');
    assert(v2_255.direct_v2 == 20 && v2_127.direct_v2 == 17 && ...
        v2_255.bch_direct_cyc == 3 && v2_127.bch_direct_cyc == 2, ...
        'v2 small-GF direct latency mismatch.');
    [ratio, kpiPass] = HardCascade.kpi_check(v2_255.direct_v2, v2_255.pure_rs);
    assert(kpiPass && abs(ratio - 100/19) < 1e-12, 'n255 v2 KPI check mismatch.');

    %% Both public experiment entry points execute in test-safe profile mode.
    profileMain = Work03Experiments.main(true, false);
    assert(isfield(profileMain, 'results') && isfield(profileMain.results, 'n255') && ...
        isfield(profileMain.results.n255, 'pure_rs') && ...
        numel(profileMain.results.n127.cascade_direct.fer) >= 1, ...
        'main(profile) did not return all configured methods.');
    profileV2 = Work03Experiments.main_v2(true, true);
    paths = HardCascade.paths();
    assert(profileV2.n255_cycles == 20 && profileV2.n127_cycles == 17, ...
        'main_v2(profile) returned incorrect direct-v2 cycles.');
    assert(exist(fullfile(paths.data_root, 'matlab_profile_v2_latency.json'), 'file') == 2 && ...
        exist(fullfile(paths.figures_root, 'matlab_profile_latency_bars_v2.png'), 'file') == 2 && ...
        exist(fullfile(paths.figures_root, 'matlab_profile_latency_evolution.pdf'), 'file') == 2, ...
        'main_v2 did not write Work 03-local assets.');

    report = HardCascade.selftest();
    assert(all([report.pass]), 'HardCascade.selftest reported a failure.');
    fprintf('WORK03_MATLAB_TESTS_OK\n');
end
