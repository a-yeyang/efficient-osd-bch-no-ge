function [c_hat, ok] = bch_decode_direct(gf, r, lut_A)
    % BCH t=2 Direct root finding (Lagendijk 2026 §III-A)
    % r: length-n binary received word
    % lut_A: precomputed LUT from build_lut_A(gf)
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

    % 2 errors: transform Λ_monic → A(Y) = Y² + Y + k
    % k = (S1³ + S3) / S1³
    numer = bitxor(S1_cubed, S3);
    k_lut = gf_div(gf, numer, S1_cubed);

    if ~lut_A.valid(k_lut + 1)
        c_hat = r; ok = false; return;
    end

    Y1 = lut_A.roots(k_lut + 1, 1);
    Y2 = lut_A.roots(k_lut + 1, 2);

    % Roots of Λ_monic: X_i = S1 * Y_i, and X_i = α^{p_i} → p_i = log(X_i)
    X1 = gf_mul(gf, S1, Y1);
    X2 = gf_mul(gf, S1, Y2);
    if X1 == 0 || X2 == 0
        c_hat = r; ok = false; return;
    end
    p1 = mod(gf.LOG(X1 + 1), n);
    p2 = mod(gf.LOG(X2 + 1), n);
    if p1 == p2
        c_hat = r; ok = false; return;
    end

    c_hat = r;
    c_hat(p1+1) = 1 - c_hat(p1+1);
    c_hat(p2+1) = 1 - c_hat(p2+1);
    ok = true;
end
