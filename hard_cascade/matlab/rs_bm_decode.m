function [c, ok] = rs_bm_decode(rs, r)
    % Berlekamp-Massey + Chien + Forney hard-decision RS decode.
    % Port of rs_bch_cascade_matlab/RSCode.m bmDecode (== cascade_src/rs_code.py).
    gf = rs.gf; n = rs.n; t = rs.t;
    r = double(r(:).');

    % 1. Syndromes S_i = r(alpha^i), i = 1..2t.
    S = zeros(1, 2*t + 1);
    nz = find(r ~= 0);
    if isempty(nz)
        c = r; ok = true; return;
    end
    j0 = nz - 1;
    r_log = gf.LOG(r(nz) + 1);
    for i = 1:(2*t)
        exps = mod(i * j0, gf.n);
        prod_log = mod(r_log + exps, gf.n);
        prod = gf.EXP(prod_log + 1);
        acc = 0;
        for q = 1:numel(prod), acc = bitxor(acc, prod(q)); end
        S(i + 1) = acc;
    end

    if ~any(S(2:end))
        c = r; ok = true; return;
    end

    % 2. Berlekamp-Massey.
    L = 0;
    Lam = 1;
    B = 1;
    b = 1;
    m_shift = 1;
    for kk = 1:(2*t)
        delta = S(kk + 1);
        for i = 1:L
            if (i+1) <= numel(Lam) && Lam(i+1) ~= 0
                delta = bitxor(delta, gf_mul(gf, Lam(i+1), S(kk - i + 1)));
            end
        end
        if delta == 0
            m_shift = m_shift + 1;
        else
            coef = gf_div(gf, delta, b);
            xmB = [zeros(1, m_shift), B];
            new_len = max(numel(Lam), numel(xmB));
            T = [Lam, zeros(1, new_len - numel(Lam))];
            xmB = [xmB, zeros(1, new_len - numel(xmB))];
            for i = 1:new_len
                T(i) = bitxor(T(i), gf_mul(gf, coef, xmB(i)));
            end
            if 2*L <= kk - 1
                L_new = kk - L;
                B = Lam;
                b = delta;
                Lam = T;
                L = L_new;
                m_shift = 1;
            else
                Lam = T;
                m_shift = m_shift + 1;
            end
        end
    end

    % 3. Chien search: roots of Lam.
    Lam_deg0 = find(Lam ~= 0) - 1;
    if isempty(Lam_deg0)
        c = r; ok = false; return;
    end
    log_lam = gf.LOG(Lam(Lam_deg0 + 1) + 1);
    err_positions = [];
    for i = 0:n-1
        exps = mod((n - i) * Lam_deg0, gf.n);
        prod_log = mod(log_lam + exps, gf.n);
        prod = gf.EXP(prod_log + 1);
        val = 0;
        for q = 1:numel(prod), val = bitxor(val, prod(q)); end
        if val == 0
            err_positions(end+1) = i; %#ok<AGROW>
        end
    end

    if numel(err_positions) ~= L || L > t
        c = r; ok = false; return;
    end

    % 4. Forney. Omega(x) = [S(x) Lam(x)] mod x^{2t+1}.
    Sx = [0, S(2:2*t+1)];
    Omega = zeros(1, 2*t + 1);
    for i = 0:numel(Sx)-1
        if Sx(i+1) == 0, continue; end
        for jj = 0:numel(Lam)-1
            if (i + jj) <= 2*t && Lam(jj+1) ~= 0
                Omega(i + jj + 1) = bitxor(Omega(i + jj + 1), ...
                    gf_mul(gf, Sx(i+1), Lam(jj+1)));
            end
        end
    end
    Lam_prime = zeros(1, numel(Lam));
    for i = 1:numel(Lam)-1
        if mod(i, 2) == 1
            Lam_prime(i) = Lam(i + 1);
        end
    end

    c = r; ok = true;
    for pidx = 1:numel(err_positions)
        p = err_positions(pidx);
        omega_val = 0;
        for i = 0:numel(Omega)-1
            if Omega(i+1) == 0, continue; end
            omega_val = bitxor(omega_val, ...
                gf_mul(gf, Omega(i+1), gf.EXP(mod(i*(n-p), gf.n) + 1)));
        end
        lam_prime_val = 0;
        for i = 0:numel(Lam_prime)-1
            if Lam_prime(i+1) == 0, continue; end
            lam_prime_val = bitxor(lam_prime_val, ...
                gf_mul(gf, Lam_prime(i+1), gf.EXP(mod(i*(n-p), gf.n) + 1)));
        end
        if lam_prime_val == 0
            c = r; ok = false; return;
        end
        err_val = gf_div(gf, omega_val, lam_prime_val);
        alpha_p = gf.EXP(mod(p, gf.n) + 1);
        err_val = gf_mul(gf, alpha_p, err_val);
        c(p + 1) = bitxor(c(p + 1), err_val);
    end
end
