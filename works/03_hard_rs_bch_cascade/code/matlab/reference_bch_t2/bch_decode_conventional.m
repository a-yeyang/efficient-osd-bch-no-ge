function [c_hat, ok] = bch_decode_conventional(gf, r)
    % BCH t=2 Conventional decoder: BM (simplified) + Chien search
    % r: length-n binary received word (row vector of 0/1)
    % Returns corrected codeword c_hat and success flag ok
    n = gf.n;
    [S1, S3] = bch_syndromes(gf, r);
    S1_cubed = gf_mul(gf, gf_mul(gf, S1, S1), S1);

    if S1 == 0 && S3 == 0
        c_hat = r; ok = true; return;
    end
    if S1 == 0 && S3 ~= 0
        c_hat = r; ok = false; return;
    end

    if S1_cubed == S3
        % 1 error: p = log(S1)
        log_s1 = gf.LOG(S1 + 1);
        p = mod(log_s1, n);
        c_hat = r;
        c_hat(p+1) = 1 - c_hat(p+1);
        ok = true; return;
    end

    % 2 errors: monic ELP = x^2 + S1 x + c0
    c0 = gf_div(gf, bitxor(S1_cubed, S3), S1);

    % Chien-style search: find i s.t. α^{2i} ⊕ S1·α^i ⊕ c0 = 0
    err_pos = [];
    for i = 0:n-1
        X = gf.EXP(i+1);
        X_sq = gf.EXP(mod(2*i, n)+1);
        val = bitxor(bitxor(X_sq, gf_mul(gf, S1, X)), c0);
        if val == 0
            err_pos = [err_pos, i];
            if length(err_pos) == 2, break; end
        end
    end
    if length(err_pos) ~= 2
        c_hat = r; ok = false; return;
    end
    c_hat = r;
    for i = 1:length(err_pos)
        p = err_pos(i);
        c_hat(p+1) = 1 - c_hat(p+1);
    end
    ok = true;
end
