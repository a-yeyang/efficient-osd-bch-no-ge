function ok = smoke_test_v3()
    % Correctness smoke test for the v3 hard-decision cascade machinery.
    % Mirrors the Python smoke tests in docs/report_v3.tex Sec. 2.3:
    %   - RS(544,514,t=15)/GF(2^10): inject <=15 symbol errors, must recover.
    %   - BCH(144,136,t=1)/GF(2^8): inject 1 bit error, must recover.
    %   - BCH(255,239,t=2)/GF(2^8): inject <=2 bit errors, must recover.
    %   - Full cascade round-trip (both configs) with zero channel noise.
    here = fileparts(mfilename('fullpath'));
    addpath(here);
    rng(1);
    ok = true;

    fprintf('=== smoke_test_v3: 正确性自检 ===\n');

    % --- RS(544,514,t=15)/GF(2^10), shortened from RS(1023,993) -----------
    srs = shortened_rs_init(10, 544, 514);
    n_pass = 0;
    for trial = 1:20
        msg = randi([0, 1023], 1, srs.k);
        c = shortened_rs_encode(srs, msg);
        r = c;
        n_err = randi([0, srs.t]);
        pos = randperm(srs.n, n_err);
        for p = pos
            wrong = mod(r(p) + 1 + randi(1022), 1023);  % any symbol != r(p)
            r(p) = wrong;
        end
        [c_hat, dok] = shortened_rs_decode(srs, r);
        if dok && all(c_hat == c)
            n_pass = n_pass + 1;
        end
    end
    fprintf('RS(544,514,t=15)/GF(2^10): %d/20 通过 (<=15 符号错)\n', n_pass);
    ok = ok && (n_pass == 20);

    % --- BCH(144,136,t=1)/GF(2^8), shortened from BCH(255,247) -------------
    bch1 = bch_init(8, 1, 144, 136);
    n_pass = 0;
    for trial = 1:20
        msg = randi([0, 1], 1, bch1.k);
        cw = bch_encode(bch1, msg);
        r = cw;
        p = randi(bch1.n);
        r(p) = 1 - r(p);
        [cw_hat_c, ok_c] = bch_decode_conventional_shortened(bch1, r);
        [cw_hat_d, ok_d] = bch_decode_direct_shortened(bch1, r);
        if ok_c && ok_d && all(cw_hat_c == cw) && all(cw_hat_d == cw)
            n_pass = n_pass + 1;
        end
    end
    fprintf('BCH(144,136,t=1)/GF(2^8): %d/20 通过 (1 bit 错, conv+direct)\n', n_pass);
    ok = ok && (n_pass == 20);

    % --- BCH(255,239,t=2)/GF(2^8), natural (unshortened) -------------------
    bch2 = bch_init(8, 2, 255, 239);
    n_pass = 0;
    for trial = 1:20
        msg = randi([0, 1], 1, bch2.k);
        cw = bch_encode(bch2, msg);
        r = cw;
        n_err = randi([0, 2]);
        pos = randperm(bch2.n, n_err);
        r(pos) = 1 - r(pos);
        [cw_hat_c, ok_c] = bch_decode_conventional_shortened(bch2, r);
        [cw_hat_d, ok_d] = bch_decode_direct_shortened(bch2, r);
        if ok_c && ok_d && all(cw_hat_c == cw) && all(cw_hat_d == cw)
            n_pass = n_pass + 1;
        end
    end
    fprintf('BCH(255,239,t=2)/GF(2^8): %d/20 通过 (<=2 bit 错, conv+direct)\n', n_pass);
    ok = ok && (n_pass == 20);

    % --- Full cascade round-trip, zero noise --------------------------------
    for cname = {'cfg1_kp4', 'cfg2_255'}
        cfg = cascade_config(cname{1});
        codec = cascade_init(cfg, 'direct');
        msg = randi([0, (2^cfg.m_rs)-1], 1, cfg.k_rs);
        coded = cascade_encode(codec, msg);
        [msg_hat, ok_rs, ~] = cascade_decode(codec, coded); %#ok<ASGLU>
        pass = ok_rs && all(msg_hat == msg);
        fprintf('级联往返 (零噪声) %s: %s\n', cfg.name, ternary(pass, 'PASS', 'FAIL'));
        ok = ok && pass;
    end

    if ok
        fprintf('=== smoke_test_v3: 全部通过 ===\n');
    else
        fprintf('=== smoke_test_v3: 存在失败项 ===\n');
    end
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
