classdef RSCode
%RSCODE  Narrow-sense primitive Reed-Solomon code over GF(2^m).
%
%   Port of cascade_src/rs_code.py (RSCode class).
%
%   Systematic RS(n,k) with n = 2^m - 1. Codeword layout matches the Python
%   reference EXACTLY:  c = [parity(1..n-k), message(1..k)]  (parity first).
%   Decoders: bmDecode (Berlekamp-Massey + Chien + Forney, hard) and
%   lccBrDecode (Chase-BM soft on eta least-reliable positions).
%
%   Symbols are 0-based field-element values 0..2^m-1, stored as row vectors.

    properties
        gf          % GF object
        m           % GF exponent
        n           % 2^m - 1
        k           % message symbols
        d           % n - k + 1
        t           % (d-1)/2
        alpha_pow   % alpha^0..alpha^{n-1}
        g_poly      % generator polynomial (little-endian, coeffs in GF(2^m))
    end

    methods
        function obj = RSCode(m, k)
            obj.gf = GF(m);
            obj.m = m;
            obj.n = obj.gf.n;
            obj.k = k;
            obj.d = obj.n - obj.k + 1;
            obj.t = floor((obj.d - 1) / 2);

            obj.alpha_pow = obj.gf.EXP(1:obj.n);  % alpha^0..alpha^{n-1}

            % g(x) = prod_{i=1..2t} (x - alpha^i), char 2 => minus = plus.
            g = 1;
            for i = 1:(2*obj.t)
                root = obj.gf.EXP(i + 1);   % alpha^i  (0-based exponent i)
                g = obj.gf.polyMul(g, [root, 1]);
            end
            obj.g_poly = g;
            assert(numel(obj.g_poly) - 1 == obj.n - obj.k, ...
                'g_poly degree mismatch');
        end

        % -----------------------------------------------------------------
        function c = encodeSystematic(obj, msg)
            % Systematic encode: c = [parity(n-k), msg(k)].
            assert(numel(msg) == obj.k);
            gf = obj.gf;
            nk = obj.n - obj.k;
            % dividend = msg * x^{n-k} : [0,...,0(n-k), msg_0,...,msg_{k-1}]
            dividend = [zeros(1, nk), double(msg(:).')];
            [~, rem] = gf.polyDivmod(dividend, obj.g_poly);
            parity = zeros(1, nk);
            parity(1:numel(rem)) = rem;
            c = zeros(1, obj.n);
            c(1:nk) = parity;
            c(nk+1:end) = msg(:).';
        end

        function msg = extractMessage(obj, c)
            msg = c(obj.n - obj.k + 1 : end);
        end

        % -----------------------------------------------------------------
        function [c, ok] = bmDecode(obj, r, counters)
            % Berlekamp-Massey + Chien + Forney hard-decision decode.
            if nargin < 3 || isempty(counters), counters = OpCounters(); end
            gf = obj.gf; n = obj.n; t = obj.t;
            r = double(r(:).');

            % 1. Syndromes S_i = r(alpha^i), i = 1..2t.
            S = zeros(1, 2*t + 1);          % S(i+1) = S_i, S(1)=S_0 unused
            nz = find(r ~= 0);              % 1-based positions of nonzero r
            if isempty(nz)
                c = r; ok = true; return;
            end
            j0 = nz - 1;                    % 0-based symbol indices
            r_log = gf.LOG(r(nz) + 1);      % logs of nonzero symbols
            for i = 1:(2*t)
                exps = mod(i * j0, gf.n);           % alpha^{i*j} exponents
                prod_log = mod(r_log + exps, gf.n); % log(r_j) + log(alpha^{ij})
                prod = gf.EXP(prod_log + 1);
                acc = 0;
                for q = 1:numel(prod), acc = bitxor(acc, prod(q)); end
                S(i + 1) = acc;
            end
            counters.f2m = counters.f2m + numel(nz) * 2 * t;

            if ~any(S(2:end))
                c = r; ok = true; return;
            end

            % 2. Berlekamp-Massey.
            L = 0;
            Lam = 1;        % little-endian (Lam(1)=coeff of x^0)
            B = 1;
            b = 1;
            m_shift = 1;
            for kk = 1:(2*t)
                delta = S(kk + 1);
                for i = 1:L
                    if (i+1) <= numel(Lam) && Lam(i+1) ~= 0
                        delta = bitxor(delta, gf.mul(Lam(i+1), S(kk - i + 1)));
                        counters.f2m = counters.f2m + 1;
                    end
                end
                if delta == 0
                    m_shift = m_shift + 1;
                else
                    coef = gf.div(delta, b);
                    counters.f2m = counters.f2m + 1;
                    xmB = [zeros(1, m_shift), B];
                    new_len = max(numel(Lam), numel(xmB));
                    T = [Lam, zeros(1, new_len - numel(Lam))];
                    xmB = [xmB, zeros(1, new_len - numel(xmB))];
                    for i = 1:new_len
                        T(i) = bitxor(T(i), gf.mul(coef, xmB(i)));
                        counters.f2m = counters.f2m + 1;
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
            Lam_deg0 = find(Lam ~= 0) - 1;   % 0-based degrees of nonzero terms
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
            counters.f2m = counters.f2m + n * numel(Lam_deg0);

            if numel(err_positions) ~= L || L > t
                c = r; ok = false; return;
            end

            % 4. Forney. Omega(x) = [S(x) Lam(x)] mod x^{2t+1}.
            Sx = [0, S(2:2*t+1)];   % Sx(1)=0, Sx(i+1)=S_i
            Omega = zeros(1, 2*t + 1);
            for i = 0:numel(Sx)-1
                if Sx(i+1) == 0, continue; end
                for jj = 0:numel(Lam)-1
                    if (i + jj) <= 2*t && Lam(jj+1) ~= 0
                        Omega(i + jj + 1) = bitxor(Omega(i + jj + 1), ...
                            gf.mul(Sx(i+1), Lam(jj+1)));
                        counters.f2m = counters.f2m + 1;
                    end
                end
            end
            % Lam'(x): odd-degree terms only (char 2).
            Lam_prime = zeros(1, numel(Lam));
            for i = 1:numel(Lam)-1     % i is 0-based degree index -> Lam(i+1)
                if mod(i, 2) == 1
                    Lam_prime(i) = Lam(i + 1);   % shifts down by one degree
                end
            end

            c = r; ok = true;
            for pidx = 1:numel(err_positions)
                p = err_positions(pidx);
                omega_val = 0;
                for i = 0:numel(Omega)-1
                    if Omega(i+1) == 0, continue; end
                    omega_val = bitxor(omega_val, ...
                        gf.mul(Omega(i+1), gf.EXP(mod(i*(n-p), gf.n) + 1)));
                    counters.f2m = counters.f2m + 1;
                end
                lam_prime_val = 0;
                for i = 0:numel(Lam_prime)-1
                    if Lam_prime(i+1) == 0, continue; end
                    lam_prime_val = bitxor(lam_prime_val, ...
                        gf.mul(Lam_prime(i+1), gf.EXP(mod(i*(n-p), gf.n) + 1)));
                    counters.f2m = counters.f2m + 1;
                end
                if lam_prime_val == 0
                    c = r; ok = false; return;
                end
                err_val = gf.div(omega_val, lam_prime_val);
                alpha_p = gf.EXP(mod(p, gf.n) + 1);
                err_val = gf.mul(alpha_p, err_val);
                counters.f2m = counters.f2m + 1;
                c(p + 1) = bitxor(c(p + 1), err_val);
            end
        end

        % -----------------------------------------------------------------
        function [best_c, ok] = lccBrDecode(obj, r_soft, reliability, eta, counters)
            % Chase-BM soft decode: flip 2^eta patterns on the eta LRPs, BM each,
            % score by summed reliability where c_hat == r_soft, keep best.
            if nargin < 5 || isempty(counters), counters = OpCounters(); end
            r_soft = double(r_soft(:).');
            reliability = double(reliability(:).');

            [~, order] = sort(reliability, 'ascend');
            lrp = order(1:eta);            % 1-based least-reliable positions

            best_c = [];
            best_score = -inf;
            n_tvs = 0;

            for mask = 0:(2^eta - 1)
                n_tvs = n_tvs + 1;
                r_test = r_soft;
                for i = 1:eta
                    if bitand(mask, bitshift(1, i-1)) ~= 0
                        r_test(lrp(i)) = bitxor(r_test(lrp(i)), 1);
                    end
                end
                [c_hat, okk] = obj.bmDecode(r_test, counters);
                if ~okk, continue; end
                match = (c_hat == r_soft);
                score = sum(reliability(match));
                if score > best_score
                    best_score = score;
                    best_c = c_hat;
                end
            end

            if isempty(best_c)
                [c_hat, okk] = obj.bmDecode(r_soft, counters);
                if okk, best_c = c_hat; ok = true;
                else, best_c = r_soft; ok = false; end
                return;
            end
            counters.n_tvs = counters.n_tvs + n_tvs;
            ok = true;
        end
    end
end
