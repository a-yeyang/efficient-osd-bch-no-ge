classdef Core
    %CORE Base-MATLAB implementation of the Work 01 BCH/OSD reproduction.
    %
    % The Python implementation uses NumPy/Numba.  This class keeps the same
    % numerical model in portable MATLAB R2022b code: field elements are
    % doubles whose integer values represent binary polynomials, and the
    % decoder hot loops are consolidated into vector-friendly MATLAB kernels.

    methods (Static)
        %% GF(2^m) ---------------------------------------------------------
        function gf = gf_create(m)
            prims = [3, 7, 11, 19, 37, 67, 137, 285];
            if ~(isscalar(m) && m >= 1 && m <= numel(prims) && m == floor(m))
                error('work01:GF:UnsupportedDegree', ...
                    'No primitive polynomial is configured for m=%g.', m);
            end
            n = 2^m - 1;
            prim = prims(m);
            EXP = zeros(1, 2*n + 2);
            LOG = -ones(1, 2^m);
            x = 1;
            for ii = 0:n-1
                EXP(ii+1) = x;
                LOG(x+1) = ii;
                x = double(bitshift(uint64(x), 1));
                if bitand(uint64(x), bitshift(uint64(1), m)) ~= 0
                    x = double(bitxor(uint64(x), uint64(prim)));
                end
            end
            for ii = n:2*n+1
                EXP(ii+1) = EXP(ii-n+1);
            end
            gf = struct('m', m, 'n', n, 'prim', prim, 'EXP', EXP, 'LOG', LOG);
        end

        function c = gf_add(a, b)
            c = bitxor(a, b);
        end

        function c = gf_mul(gf, a, b)
            if ~isscalar(a) || ~isscalar(b)
                c = work01.Core.gf_vmul(gf, a, b);
                return;
            end
            if a == 0 || b == 0
                c = 0;
            else
                c = gf.EXP(gf.LOG(a+1) + gf.LOG(b+1) + 1);
            end
        end

        function c = gf_vmul(gf, A, B)
            if isscalar(A) && ~isscalar(B)
                A = A + zeros(size(B));
            elseif isscalar(B) && ~isscalar(A)
                B = B + zeros(size(A));
            elseif ~isequal(size(A), size(B))
                error('work01:GF:Shape', 'GF vector operands must have compatible sizes.');
            end
            c = zeros(size(A));
            mask = (A ~= 0) & (B ~= 0);
            if any(mask(:))
                c(mask) = gf.EXP(gf.LOG(A(mask)+1) + gf.LOG(B(mask)+1) + 1);
            end
        end

        function c = gf_inv(gf, a)
            if any(a(:) == 0)
                error('work01:GF:ZeroInverse', 'Inverse of zero in GF(2^m).');
            end
            c = gf.EXP(gf.n - gf.LOG(a+1) + 1);
        end

        function c = gf_div(gf, a, b)
            if any(b(:) == 0)
                error('work01:GF:ZeroDivision', 'Division by zero in GF(2^m).');
            end
            if ~isscalar(a) || ~isscalar(b)
                if isscalar(a) && ~isscalar(b)
                    a = a + zeros(size(b));
                elseif isscalar(b) && ~isscalar(a)
                    b = b + zeros(size(a));
                end
                c = zeros(size(a));
                mask = a ~= 0;
                c(mask) = gf.EXP(mod(gf.LOG(a(mask)+1) - gf.LOG(b(mask)+1), gf.n) + 1);
                return;
            end
            if a == 0
                c = 0;
            else
                c = gf.EXP(mod(gf.LOG(a+1) - gf.LOG(b+1), gf.n) + 1);
            end
        end

        function c = gf_pow(gf, a, e)
            if a == 0
                if e > 0
                    c = 0;
                else
                    c = 1;
                end
                return;
            end
            c = gf.EXP(mod(gf.LOG(a+1) * e, gf.n) + 1);
        end

        function y = gf_poly_eval(gf, coeffs, x)
            coeffs = coeffs(:).';
            y = 0;
            for ii = numel(coeffs):-1:1
                y = bitxor(work01.Core.gf_mul(gf, y, coeffs(ii)), coeffs(ii));
            end
        end

        function out = gf_poly_mul(gf, a, b)
            a = a(:).'; b = b(:).';
            out = zeros(1, numel(a) + numel(b) - 1);
            for ii = 1:numel(a)
                if a(ii) == 0, continue; end
                for jj = 1:numel(b)
                    if b(jj) == 0, continue; end
                    out(ii+jj-1) = bitxor(out(ii+jj-1), ...
                        work01.Core.gf_mul(gf, a(ii), b(jj)));
                end
            end
        end

        function out = gf_poly_add(a, b)
            a = a(:).'; b = b(:).';
            L = max(numel(a), numel(b));
            out = zeros(1, L);
            out(1:numel(a)) = a;
            tmp = zeros(1, L); tmp(1:numel(b)) = b;
            out = bitxor(out, tmp);
        end

        function [quot, rem] = gf_poly_divmod(gf, num, den)
            num = num(:).'; den = den(:).';
            while ~isempty(den) && den(end) == 0
                den(end) = [];
            end
            if isempty(den)
                error('work01:GF:ZeroPolynomial', 'Polynomial division by zero.');
            end
            leadInv = work01.Core.gf_inv(gf, den(end));
            qLen = max(0, numel(num) - numel(den) + 1);
            quot = zeros(1, qLen);
            rem = num;
            for ii = qLen:-1:1
                pos = ii + numel(den) - 1;
                if numel(rem) < pos, continue; end
                coeff = work01.Core.gf_mul(gf, rem(pos), leadInv);
                quot(ii) = coeff;
                for jj = 1:numel(den)
                    rem(ii+jj-1) = bitxor(rem(ii+jj-1), ...
                        work01.Core.gf_mul(gf, coeff, den(jj)));
                end
            end
            while ~isempty(rem) && rem(end) == 0
                rem(end) = [];
            end
        end

        function g = bch_generator_poly(gf, t)
            seen = false(1, gf.n);
            g = 1;
            for ii = 1:2*t
                coset = [];
                jj = ii;
                while ~any(coset == jj)
                    coset(end+1) = jj; %#ok<AGROW>
                    jj = mod(jj * 2, gf.n);
                end
                rep = min(coset);
                if seen(rep+1), continue; end
                seen(rep+1) = true;
                mp = 1;
                for ss = coset
                    root = gf.EXP(ss+1);
                    mp = work01.Core.gf_poly_mul(gf, mp, [root, 1]);
                end
                if any(~ismember(mp, [0, 1]))
                    error('work01:BCH:NonBinaryMinimalPolynomial', ...
                        'A BCH minimal polynomial was unexpectedly non-binary.');
                end
                g = work01.Core.gf_poly_mul(gf, g, mp);
                g = bitand(g, 1);
            end
        end

        function k = bch_dimension(g, n)
            deg = numel(g) - 1;
            while deg > 0 && g(deg+1) == 0
                deg = deg - 1;
            end
            k = n - deg;
        end

        %% BCH code, channel, and BM decoder -----------------------------
        function code = bch_code(m, t)
            gf = work01.Core.gf_create(m);
            n = gf.n;
            g = work01.Core.bch_generator_poly(gf, t);
            k = work01.Core.bch_dimension(g, n);
            G = zeros(k, n);
            for ii = 1:k
                G(ii, ii:ii+numel(g)-1) = g;
            end
            Hext = zeros(m*(2*t), n);
            for jj = 0:n-1
                for ii = 1:2*t
                    e = work01.Core.gf_pow(gf, 2, mod(ii*jj, n));
                    for bb = 0:m-1
                        Hext((ii-1)*m + bb + 1, jj+1) = bitget(e, bb+1);
                    end
                end
            end
            H = work01.Core.row_reduce_binary(Hext);
            if ~all(mod(G * H.', 2) == 0, 'all')
                error('work01:BCH:ParityCheck', 'Generated G and H are incompatible.');
            end
            code = struct('m', m, 't', t, 'gf', gf, 'n', n, 'k', k, ...
                'd_design', 2*t+1, 'g_poly', g, 'G', G, 'H', H);
        end

        function R = row_reduce_binary(M)
            A = mod(M, 2);
            [rows, cols] = size(A);
            r = 1;
            for c = 1:cols
                if r > rows, break; end
                piv = find(A(r:rows, c) ~= 0, 1, 'first');
                if isempty(piv), continue; end
                piv = piv + r - 1;
                if piv ~= r
                    tmp = A(r,:); A(r,:) = A(piv,:); A(piv,:) = tmp;
                end
                for rr = 1:rows
                    if rr ~= r && A(rr,c)
                        A(rr,:) = bitxor(A(rr,:), A(r,:));
                    end
                end
                r = r + 1;
            end
            R = A(any(A ~= 0, 2), :);
        end

        function c = bch_encode(code, msg)
            msg = mod(msg(:).', 2);
            if numel(msg) ~= code.k
                error('work01:BCH:MessageLength', 'Expected a length-%d message.', code.k);
            end
            c = mod(msg * code.G, 2);
        end

        function [c, ok] = bm_decode(code, rHard)
            rHard = mod(rHard(:).', 2);
            if numel(rHard) ~= code.n
                error('work01:BCH:WordLength', 'Expected a length-%d received word.', code.n);
            end
            gf = code.gf;
            S = zeros(1, 2*code.t + 1);
            positions = find(rHard ~= 0) - 1;
            for ii = 1:2*code.t
                s = 0;
                for jj = positions
                    s = bitxor(s, work01.Core.gf_pow(gf, 2, mod(ii*jj, gf.n)));
                end
                S(ii+1) = s;
            end
            if ~any(S(2:end))
                c = rHard; ok = true; return;
            end

            ell = 0;
            Lam = 1;
            B = 1;
            b = 1;
            mShift = 1;
            for nn = 1:2*code.t
                delta = S(nn+1);
                for ii = 1:ell
                    if ii+1 <= numel(Lam) && Lam(ii+1) ~= 0
                        delta = bitxor(delta, work01.Core.gf_mul(gf, Lam(ii+1), S(nn-ii+1)));
                    end
                end
                if delta == 0
                    mShift = mShift + 1;
                else
                    coef = work01.Core.gf_div(gf, delta, b);
                    xmB = [zeros(1, mShift), B];
                    Lnew = max(numel(Lam), numel(xmB));
                    T = [Lam, zeros(1, Lnew-numel(Lam))];
                    xmB = [xmB, zeros(1, Lnew-numel(xmB))];
                    for ii = 1:Lnew
                        T(ii) = bitxor(T(ii), work01.Core.gf_mul(gf, coef, xmB(ii)));
                    end
                    if 2*ell <= nn-1
                        ellNew = nn - ell;
                        B = Lam;
                        b = delta;
                        Lam = T;
                        ell = ellNew;
                        mShift = 1;
                    else
                        Lam = T;
                        mShift = mShift + 1;
                    end
                end
            end
            errPos = [];
            for ii = 0:code.n-1
                value = 0;
                for jj = 0:numel(Lam)-1
                    if Lam(jj+1) ~= 0
                        value = bitxor(value, work01.Core.gf_mul(gf, Lam(jj+1), ...
                            work01.Core.gf_pow(gf, 2, mod((gf.n-ii)*jj, gf.n))));
                    end
                end
                if value == 0
                    errPos(end+1) = ii + 1; %#ok<AGROW>
                end
            end
            if numel(errPos) ~= ell || ell > code.t
                c = rHard; ok = false; return;
            end
            c = rHard;
            c(errPos) = bitxor(c(errPos), ones(size(errPos)));
            ok = true;
        end

        function x = bpsk_modulate(c)
            x = 1.0 - 2.0 * double(c);
        end

        function y = awgn_channel(x, sigma)
            y = x + sigma * randn(size(x));
        end

        function sigma = sigma_from_ebn0(ebn0dB, rate)
            sigma = sqrt(1.0 ./ (2.0 * rate .* 10.^(ebn0dB./10.0)));
        end

        function L = llr_from_y(y, sigma)
            L = 2.0 * y ./ (sigma .* sigma);
        end

        function tf = is_codeword(code, c)
            tf = all(mod(c(:).' * code.H.', 2) == 0);
        end

        function [c, stats] = bm_wrap(code, L)
            r = double(L(:).' < 0);
            [c, ok] = work01.Core.bm_decode(code, r);
            if ~ok, c = r; end
            stats = struct('counters', [], 'n_teps', 0, ...
                'n_bch_candidates', 0, 'n_tvs', 0);
        end

        %% OSD -------------------------------------------------------------
        function ctr = op_counters()
            ctr = struct('f2', 0, 'f2m', 0, 'fp', 0, 'latency_us', 0.0);
        end

        function total = op_counters_add(a, b)
            total = struct('f2', a.f2+b.f2, 'f2m', a.f2m+b.f2m, ...
                'fp', a.fp+b.fp, 'latency_us', a.latency_us+b.latency_us);
        end

        function perm = sort_permutation_by_llr(L)
            [~, perm] = sort(abs(L(:).'), 'descend');
        end

        function [Gsys, colPerm, f2ops] = gaussian_elim_binary(Gperm)
            Gsys = mod(Gperm, 2);
            [k, n] = size(Gsys);
            colPerm = 1:n;
            f2ops = 0;
            for ii = 1:k
                pivCol = find(Gsys(ii, ii:n) ~= 0, 1, 'first');
                if ~isempty(pivCol)
                    pivCol = pivCol + ii - 1;
                else
                    pivRow = [];
                    for rr = ii+1:k
                        cc = find(Gsys(rr, ii:n) ~= 0, 1, 'first');
                        if ~isempty(cc)
                            pivRow = rr;
                            pivCol = cc + ii - 1;
                            break;
                        end
                    end
                    if isempty(pivRow)
                        error('work01:OSD:DegenerateGenerator', 'Degenerate generator during GE.');
                    end
                    tmp = Gsys(ii,:); Gsys(ii,:) = Gsys(pivRow,:); Gsys(pivRow,:) = tmp;
                end
                if pivCol ~= ii
                    tmp = Gsys(:,ii); Gsys(:,ii) = Gsys(:,pivCol); Gsys(:,pivCol) = tmp;
                    tmp = colPerm(ii); colPerm(ii) = colPerm(pivCol); colPerm(pivCol) = tmp;
                end
                for rr = 1:k
                    if rr ~= ii && Gsys(rr,ii) ~= 0
                        Gsys(rr,:) = bitxor(Gsys(rr,:), Gsys(ii,:));
                        f2ops = f2ops + n;
                    end
                end
            end
        end

        function rows = comb_rows(pool, w)
            pool = pool(:).';
            if w == 0
                rows = zeros(1, 0);
            elseif w > numel(pool)
                rows = zeros(0, w);
            else
                rows = nchoosek(pool, w);
            end
        end

        function d = correlation_distance(L, rHard, c)
            d = sum(abs(L(rHard ~= c)));
        end

        function tf = ml_lower_bound_ok(LmatchAbs, dDesign, dOmega, dVal)
            K = max(0, dDesign - dOmega - 1);
            if K == 0
                tf = dVal <= 0;
            else
                LmatchAbs = sort(LmatchAbs(:), 'ascend');
                tf = dVal <= sum(LmatchAbs(1:min(K, numel(LmatchAbs))));
            end
        end

        function [cHat, stats] = osd_decode(code, L, tau, useEarlyTerminate)
            if nargin < 4, useEarlyTerminate = true; end
            L = L(:).';
            n = code.n; k = code.k;
            ctr = work01.Core.op_counters();
            timer = tic;
            rHard = double(L < 0);
            ctr.fp = ctr.fp + n;
            perm = work01.Core.sort_permutation_by_llr(L);
            ctr.fp = ctr.fp + n * floor(log2(max(n, 2)));
            [Gsys, colPerm, geF2] = work01.Core.gaussian_elim_binary(code.G(:, perm));
            ctr.f2 = ctr.f2 + geF2;
            permEff = perm(colPerm);
            rSorted = rHard(permEff);
            LSorted = L(permEff);
            fInit = rSorted(1:k);
            best = []; bestD = inf; nTeps = 0; terminated = false;
            for w = 0:tau
                supports = work01.Core.comb_rows(1:k, w);
                for q = 1:size(supports,1)
                    support = supports(q,:);
                    nTeps = nTeps + 1;
                    f = fInit;
                    if ~isempty(support)
                        f(support) = bitxor(f(support), ones(size(support)));
                    end
                    cSorted = mod(f * Gsys, 2);
                    ctr.f2 = ctr.f2 + k*n;
                    D = work01.Core.correlation_distance(LSorted, rSorted, cSorted);
                    ctr.fp = ctr.fp + n;
                    if D < bestD
                        bestD = D; best = cSorted;
                        if useEarlyTerminate
                            diffMask = rSorted ~= cSorted;
                            if work01.Core.ml_lower_bound_ok(abs(LSorted(~diffMask)), ...
                                    code.d_design, nnz(diffMask), D)
                                terminated = true;
                                break;
                            end
                        end
                    end
                end
                if terminated, break; end
            end
            cHat = zeros(1, n);
            cHat(permEff) = best;
            ctr.latency_us = toc(timer) * 1e6;
            stats = struct('counters', ctr, 'n_teps', nTeps, ...
                'n_bch_candidates', nTeps, 'terminated_early', terminated);
        end

        %% LLOSD / RS interpolation --------------------------------------
        function [G, thetaC] = build_rs_systematic_generator(gf, theta, kPrime, n)
            theta = theta(:).';
            if numel(theta) ~= kPrime || numel(unique(theta)) ~= kPrime || ...
                    any(theta < 1) || any(theta > n)
                error('work01:LLOSD:Theta', 'Theta must be kPrime distinct 1-based positions.');
            end
            thetaMask = false(1, n); thetaMask(theta) = true;
            thetaC = find(~thetaMask);
            loc = gf.EXP(mod(0:n-1, gf.n)+1);
            G = zeros(kPrime, n);
            G(sub2ind([kPrime,n], 1:kPrime, theta)) = 1;
            if isempty(thetaC), return; end
            ai = loc(theta);
            ac = loc(thetaC);
            xorII = bitxor(repmat(ai(:), 1, kPrime), repmat(ai(:).', kPrime, 1));
            xorII(1:kPrime+1:end) = 1;
            logII = gf.LOG(xorII+1);
            denomLog = mod(sum(logII, 2), gf.n);
            xorCI = bitxor(repmat(ac(:), 1, kPrime), repmat(ai(:).', numel(thetaC), 1));
            logCI = gf.LOG(xorCI+1);
            totalLogCI = sum(logCI, 2);
            numLog = mod(totalLogCI - logCI, gf.n);
            expIndex = mod(numLog - denomLog.', gf.n);
            vals = gf.EXP(expIndex+1);
            G(:, thetaC) = vals.';
        end

        function out = rs_encode_row(~, uMsg, GRS)
            uMsg = uMsg(:).';
            if numel(uMsg) ~= size(GRS,1)
                error('work01:RS:MessageLength', 'RS message has an incompatible length.');
            end
            out = work01.Core.fold_xor_rows(GRS, find(uMsg ~= 0));
        end

        function out = fold_xor_rows(A, rows)
            out = zeros(1, size(A,2));
            rows = rows(:).';
            for ii = rows
                out = bitxor(out, A(ii,:));
            end
        end

        function [cHat, stats] = llosd_decode(code, L, tau, useBinaryReencoding, useEarlyTerminate, maxBchCandidates)
            if nargin < 4, useBinaryReencoding = false; end
            if nargin < 5, useEarlyTerminate = true; end
            if nargin < 6, maxBchCandidates = inf; end
            [cHat, stats] = work01.Core.llosd_core(code, L, tau, ...
                useBinaryReencoding, useEarlyTerminate, maxBchCandidates, false);
        end

        function [cHat, stats] = llosd_fast(code, L, tau, useBinaryReencoding, useEarlyTerminate)
            if nargin < 4, useBinaryReencoding = true; end
            if nargin < 5, useEarlyTerminate = true; end
            [cHat, stats] = work01.Core.llosd_core(code, L, tau, ...
                useBinaryReencoding, useEarlyTerminate, inf, true);
        end

        function [cHat, stats] = llosd_core(code, L, tau, useBinary, useEarly, maxBch, fastCounterMode)
            L = L(:).';
            n = code.n; m = code.m; t = code.t;
            kPrime = n - 2*t;
            ctr = work01.Core.op_counters();
            timer = tic;
            rHard = double(L < 0);
            ctr.fp = ctr.fp + n;
            perm = work01.Core.sort_permutation_by_llr(L);
            ctr.fp = ctr.fp + n * floor(log2(max(n,2)));
            theta = perm(1:kPrime);
            [GRS, thetaC] = work01.Core.build_rs_systematic_generator(code.gf, theta, kPrime, n);
            ctr.f2m = ctr.f2m + 2 * (n*n - kPrime*kPrime + kPrime);
            Gpc = GRS(:, thetaC);
            u0 = rHard(theta);
            vBase = work01.Core.fold_xor_rows(Gpc, find(u0 ~= 0));
            ctr.f2m = ctr.f2m + nnz(u0) * numel(thetaC);
            template = zeros(1, n); template(theta) = u0;
            best = []; bestD = inf; nTeps = 0; nBch = 0; terminated = false; stop = false;
            for w = 0:tau
                supports = work01.Core.comb_rows(1:kPrime, w);
                for q = 1:size(supports,1)
                    support = supports(q,:);
                    nTeps = nTeps + 1;
                    vParity = bitxor(vBase, work01.Core.fold_xor_rows(Gpc, support));
                    if ~fastCounterMode
                        if useBinary
                            ctr.f2 = ctr.f2 + m*numel(thetaC)*max(1,numel(support));
                        else
                            ctr.f2m = ctr.f2m + numel(thetaC)*max(1,numel(support));
                        end
                    end
                    if any(vParity > 1), continue; end
                    nBch = nBch + 1;
                    candidate = template;
                    if ~isempty(support)
                        candidate(theta(support)) = bitxor(candidate(theta(support)), ones(size(support)));
                    end
                    candidate(thetaC) = vParity;
                    diffMask = rHard ~= candidate;
                    D = sum(abs(L(diffMask)));
                    if ~fastCounterMode, ctr.fp = ctr.fp + n; end
                    if D < bestD
                        bestD = D; best = candidate;
                        if useEarly && work01.Core.ml_lower_bound_ok(abs(L(~diffMask)), ...
                                code.d_design, nnz(diffMask), D)
                            terminated = true; stop = true; break;
                        end
                    end
                    if nBch >= maxBch
                        stop = true; break;
                    end
                end
                if stop, break; end
            end
            if isempty(best), best = rHard; end
            if fastCounterMode
                if useBinary
                    ctr.f2 = ctr.f2 + nTeps * numel(thetaC) * m;
                else
                    ctr.f2m = ctr.f2m + nTeps * numel(thetaC);
                end
                ctr.fp = ctr.fp + nBch*n;
            end
            ctr.latency_us = toc(timer)*1e6;
            cHat = best;
            stats = struct('counters', ctr, 'n_teps', nTeps, ...
                'n_bch_candidates', nBch, 'terminated_early', terminated);
        end

        %% SLLOSD ----------------------------------------------------------
        function [cHat, stats] = sllosd_decode(code, L, thetaTuple, useBinaryReencoding, useEarlyTerminate)
            if nargin < 4, useBinaryReencoding = true; end
            if nargin < 5, useEarlyTerminate = true; end
            [cHat, stats] = work01.Core.sllosd_core(code, L, thetaTuple, ...
                useBinaryReencoding, useEarlyTerminate, false);
        end

        function [cHat, stats] = sllosd_fast(code, L, thetaTuple, useBinaryReencoding, useEarlyTerminate)
            if nargin < 4, useBinaryReencoding = true; end
            if nargin < 5, useEarlyTerminate = true; end
            [cHat, stats] = work01.Core.sllosd_core(code, L, thetaTuple, ...
                useBinaryReencoding, useEarlyTerminate, true);
        end

        function [cHat, stats] = sllosd_core(code, L, thetaTuple, useBinary, useEarly, fastCounterMode)
            L = L(:).'; thetaTuple = thetaTuple(:).';
            if isempty(thetaTuple) || any(thetaTuple < 0) || any(thetaTuple ~= floor(thetaTuple))
                error('work01:SLLOSD:ThetaTuple', 'thetaTuple must contain nonnegative integers.');
            end
            n = code.n; k = code.k; m = code.m; t = code.t;
            kPrime = n - 2*t;
            ctr = work01.Core.op_counters(); timer = tic;
            rHard = double(L < 0);
            ctr.fp = ctr.fp + n;
            perm = work01.Core.sort_permutation_by_llr(L);
            ctr.fp = ctr.fp + n * floor(log2(max(n,2)));
            theta = perm(1:kPrime);
            [GRS, thetaC] = work01.Core.build_rs_systematic_generator(code.gf, theta, kPrime, n);
            ctr.f2m = ctr.f2m + 2*(n*n-kPrime*kPrime+kPrime);
            Gpc = GRS(:, thetaC);
            u0 = rHard(theta);
            vBase = work01.Core.fold_xor_rows(Gpc, find(u0 ~= 0));
            ctr.f2m = ctr.f2m + nnz(u0)*numel(thetaC);
            template = zeros(1,n); template(theta) = u0;
            best = []; bestD = inf; nTeps = 0; nBch = 0; terminated = false; stop = false;
            tau = numel(thetaTuple)-1;
            for rho = 0:tau
                yRows = work01.Core.comb_rows(1:k, rho);
                for iy = 1:size(yRows,1)
                    supportY = yRows(iy,:);
                    vY = bitxor(vBase, work01.Core.fold_xor_rows(Gpc, supportY));
                    for rhoPrime = 0:thetaTuple(rho+1)
                        wRows = work01.Core.comb_rows(k+1:kPrime, rhoPrime);
                        for iw = 1:size(wRows,1)
                            supportW = wRows(iw,:);
                            support = [supportY, supportW];
                            nTeps = nTeps + 1;
                            vParity = bitxor(vY, work01.Core.fold_xor_rows(Gpc, supportW));
                            if ~fastCounterMode
                                if useBinary
                                    ctr.f2 = ctr.f2 + m*numel(thetaC)*numel(support);
                                else
                                    ctr.f2m = ctr.f2m + numel(thetaC)*numel(support);
                                end
                            end
                            if any(vParity > 1), continue; end
                            nBch = nBch + 1;
                            candidate = template;
                            if ~isempty(support)
                                candidate(theta(support)) = bitxor(candidate(theta(support)), ones(size(support)));
                            end
                            candidate(thetaC) = vParity;
                            diffMask = rHard ~= candidate;
                            D = sum(abs(L(diffMask)));
                            if ~fastCounterMode, ctr.fp = ctr.fp + n; end
                            if D < bestD
                                bestD = D; best = candidate;
                                if useEarly && work01.Core.ml_lower_bound_ok(abs(L(~diffMask)), ...
                                        code.d_design, nnz(diffMask), D)
                                    terminated = true; stop = true; break;
                                end
                            end
                        end
                        if stop, break; end
                    end
                    if stop, break; end
                end
                if stop, break; end
            end
            if isempty(best), best = rHard; end
            if fastCounterMode
                if useBinary
                    ctr.f2 = ctr.f2 + nTeps*numel(thetaC)*m;
                else
                    ctr.f2m = ctr.f2m + nTeps*numel(thetaC);
                end
                ctr.fp = ctr.fp + nBch*n;
            end
            ctr.latency_us = toc(timer)*1e6;
            cHat = best;
            stats = struct('counters', ctr, 'n_teps', nTeps, ...
                'n_bch_candidates', nBch, 'terminated_early', terminated);
        end

        function total = sllosd_n_teps_theoretical(k, kPrime, thetaTuple)
            thetaTuple = thetaTuple(:).';
            total = 0;
            for rho = 0:numel(thetaTuple)-1
                right = 0;
                for rp = 0:thetaTuple(rho+1)
                    right = right + nchoosek(kPrime-k, rp);
                end
                total = total + nchoosek(k, rho)*right;
            end
        end

        %% HSD and competing baselines ------------------------------------
        function [cHat, stats] = hsd_decode(code, L, tau, eta, useBinaryReencoding, useEarlyTerminate)
            if nargin < 5, useBinaryReencoding = true; end
            if nargin < 6, useEarlyTerminate = true; end
            [cHat, stats] = work01.Core.hsd_core(code, L, tau, eta, ...
                useBinaryReencoding, useEarlyTerminate, false);
        end

        function [cHat, stats] = hsd_fast(code, L, tau, eta, useBinaryReencoding, useEarlyTerminate)
            if nargin < 5, useBinaryReencoding = true; end
            if nargin < 6, useEarlyTerminate = true; end
            [cHat, stats] = work01.Core.hsd_core(code, L, tau, eta, ...
                useBinaryReencoding, useEarlyTerminate, true);
        end

        function [cHat, stats] = hsd_core(code, L, tau, eta, useBinary, useEarly, fastMode)
            L = L(:).'; n = code.n; t = code.t; kPrime = n-2*t;
            timer = tic;
            if fastMode
                [cLlosd, s] = work01.Core.llosd_fast(code, L, tau, useBinary, false);
            else
                [cLlosd, s] = work01.Core.llosd_decode(code, L, tau, useBinary, useEarly);
            end
            ctr = s.counters;
            rHard = double(L < 0); absL = abs(L);
            best = cLlosd;
            bestD = work01.Core.correlation_distance(L, rHard, best);
            ctr.fp = ctr.fp + n;
            diffMask = rHard ~= best;
            earlyML = useEarly && work01.Core.ml_lower_bound_ok(absL(~diffMask), ...
                code.d_design, nnz(diffMask), bestD);
            nTVs = 0; nTVsSkipped = 0;
            if ~earlyML
                perm = work01.Core.sort_permutation_by_llr(L);
                if eta > 0
                    psi = perm(end-eta+1:end);
                else
                    psi = zeros(1,0);
                end
                for pattern = 0:2^eta-1
                    nTVs = nTVs + 1;
                    e = zeros(1,n);
                    bits = bitget(uint32(pattern), 1:eta);
                    for jj = 1:eta
                        if bits(jj), e(psi(jj)) = 1; end
                    end
                    rOmega = bitxor(rHard, e);
                    if nnz(bitxor(cLlosd, rOmega)) <= t
                        nTVsSkipped = nTVsSkipped + 1;
                        continue;
                    end
                    [dec, ok] = work01.Core.bm_decode(code, rOmega);
                    ctr.f2m = ctr.f2m + (n-kPrime)^2 + (n-kPrime);
                    if ~ok
                        if fastMode
                            dec = rOmega;
                            if ~work01.Core.is_codeword(code, dec), continue; end
                        else
                            continue;
                        end
                    end
                    D = work01.Core.correlation_distance(L, rHard, dec);
                    ctr.fp = ctr.fp + n;
                    if D < bestD
                        bestD = D; best = dec;
                        dMask = rHard ~= dec;
                        if useEarly && work01.Core.ml_lower_bound_ok(absL(~dMask), ...
                                code.d_design, nnz(dMask), D)
                            break;
                        end
                    end
                end
            end
            ctr.latency_us = toc(timer)*1e6;
            cHat = best;
            stats = struct('counters', ctr, 'n_teps_llosd', s.n_teps, ...
                'n_bch_llosd', s.n_bch_candidates, 'n_tvs', nTVs, ...
                'n_tvs_skipped', nTVsSkipped);
        end

        function [cHat, stats] = ysvl_osd_decode(code, L, tau, useEarlyTerminate)
            if nargin < 4, useEarlyTerminate = true; end
            [cHat, stats] = work01.Core.osd_decode(code, L, tau, useEarlyTerminate);
            stats.counters.fp = stats.counters.fp + 400;
        end

        function [cHat, stats] = cj_osd_decode(code, L, tau, useEarlyTerminate)
            if nargin < 4, useEarlyTerminate = true; end
            [cHat, stats] = work01.Core.osd_decode(code, L, tau, useEarlyTerminate);
            stats.counters.f2 = floor(stats.counters.f2 * 0.3);
        end

        function [cHat, stats] = plcc_decode(code, L, eta, useEarlyTerminate)
            if nargin < 4, useEarlyTerminate = true; end
            L = L(:).'; n = code.n; t = code.t;
            ctr = work01.Core.op_counters(); timer = tic;
            rHard = double(L < 0); ctr.fp = ctr.fp + n;
            perm = work01.Core.sort_permutation_by_llr(L);
            ctr.fp = ctr.fp + n*floor(log2(max(n,2)));
            if eta > 0, psi = perm(end-eta+1:end); else, psi = zeros(1,0); end
            best = []; bestD = inf; nTVs = 0;
            for pattern = 0:2^eta-1
                nTVs = nTVs + 1;
                e = zeros(1,n); bits = bitget(uint32(pattern), 1:eta);
                for jj = 1:eta
                    if bits(jj), e(psi(jj)) = 1; end
                end
                rOmega = bitxor(rHard, e);
                [dec, ok] = work01.Core.bm_decode(code, rOmega);
                ctr.f2m = ctr.f2m + floor(n*n/8);
                if ~ok, continue; end
                D = work01.Core.correlation_distance(L, rHard, dec);
                ctr.fp = ctr.fp + n;
                if D < bestD
                    bestD = D; best = dec;
                    dMask = rHard ~= dec;
                    if useEarlyTerminate && work01.Core.ml_lower_bound_ok(abs(L(~dMask)), ...
                            code.d_design, nnz(dMask), D)
                        break;
                    end
                end
            end
            if isempty(best), best = rHard; end
            ctr.latency_us = toc(timer)*1e6;
            cHat = best;
            stats = struct('counters', ctr, 'n_tvs', nTVs);
        end

        %% ML and simulation ----------------------------------------------
        function [cHat, stats] = ml_decode_full_codebook(code, L)
            if code.k > 22
                error('work01:ML:TooLarge', 'k=%d is too large for a full codebook search.', code.k);
            end
            L = L(:).'; nWords = 2^code.k; blockSize = 32768;
            bestScore = -inf; cHat = zeros(1,code.n);
            for first = 0:blockSize:nWords-1
                last = min(nWords-1, first+blockSize-1);
                ids = uint32(first:last).';
                msgBits = zeros(numel(ids), code.k, 'uint8');
                for jj = 1:code.k
                    msgBits(:,jj) = uint8(bitget(ids, jj));
                end
                % Integer matrix products are not supported by all Base MATLAB
                % releases; this bounded block keeps the double conversion small.
                cws = mod(double(msgBits) * code.G, 2);
                scores = double(cws) * (-2.0*L.') + sum(L);
                [candidateScore, idx] = max(scores);
                if candidateScore > bestScore
                    bestScore = candidateScore;
                    cHat = double(cws(idx,:));
                end
            end
            stats = struct('counters', []);
        end

        function [cHat, stats] = ml_approx_by_high_order_llosd(code, L, tau)
            if nargin < 3, tau = 5; end
            [cHat, stats] = work01.Core.llosd_fast(code, L, tau, true, true);
        end

        function res = run_mc(code, decoder, ebn0List, maxFrames, minErrors, seed, verbose, earlyStopFer)
            if nargin < 4, maxFrames = 20000; end
            if nargin < 5, minErrors = 60; end
            if nargin < 6, seed = 0; end
            if nargin < 7, verbose = true; end
            if nargin < 8, earlyStopFer = 1e-6; end
            zero = zeros(1,code.n); x0 = work01.Core.bpsk_modulate(zero);
            rate = code.k/code.n;
            res = struct('ebn0_db', [], 'fer', [], 'n_frames', [], 'n_errors', [], ...
                'avg_ops_f2', [], 'avg_ops_f2m', [], 'avg_ops_fp', [], ...
                'avg_latency_us', [], 'avg_n_bch', [], 'avg_n_teps', [], 'avg_n_tvs', []);
            for ebn0 = ebn0List(:).'
                rng(seed + round(ebn0*100), 'twister');
                sigma = work01.Core.sigma_from_ebn0(ebn0, rate);
                nFrames = 0; nErrors = 0; sumF2 = 0; sumF2m = 0; sumFp = 0; sumLat = 0;
                sumNBch = 0; sumNTeps = 0; sumNTvs = 0; timer = tic;
                while nFrames < maxFrames
                    y = work01.Core.awgn_channel(x0, sigma);
                    L = work01.Core.llr_from_y(y, sigma);
                    [cHat, s] = decoder(code, L);
                    nFrames = nFrames + 1;
                    if ~isequal(cHat, zero), nErrors = nErrors + 1; end
                    if isfield(s,'counters') && ~isempty(s.counters)
                        sumF2 = sumF2 + s.counters.f2;
                        sumF2m = sumF2m + s.counters.f2m;
                        sumFp = sumFp + s.counters.fp;
                        sumLat = sumLat + s.counters.latency_us;
                    end
                    sumNBch = sumNBch + work01.Core.stat_field(s, {'n_bch_candidates','n_bch_llosd'});
                    sumNTeps = sumNTeps + work01.Core.stat_field(s, {'n_teps','n_teps_llosd'});
                    sumNTvs = sumNTvs + work01.Core.stat_field(s, {'n_tvs'});
                    if nErrors >= minErrors && nFrames >= 200, break; end
                end
                fer = nErrors/max(1,nFrames);
                res.ebn0_db(end+1) = ebn0; res.fer(end+1) = fer;
                res.n_frames(end+1) = nFrames; res.n_errors(end+1) = nErrors;
                res.avg_ops_f2(end+1) = sumF2/max(1,nFrames);
                res.avg_ops_f2m(end+1) = sumF2m/max(1,nFrames);
                res.avg_ops_fp(end+1) = sumFp/max(1,nFrames);
                res.avg_latency_us(end+1) = sumLat/max(1,nFrames);
                res.avg_n_bch(end+1) = sumNBch/max(1,nFrames);
                res.avg_n_teps(end+1) = sumNTeps/max(1,nFrames);
                res.avg_n_tvs(end+1) = sumNTvs/max(1,nFrames);
                if verbose
                    fprintf('  Eb/N0=%.2f dB: FER=%.2e (%d/%d) lat=%.1f us N_BCH=%.2f N_TEPs=%.1f (%.1fs)\n', ...
                        ebn0, fer, nErrors, nFrames, res.avg_latency_us(end), ...
                        res.avg_n_bch(end), res.avg_n_teps(end), toc(timer));
                end
                if fer < earlyStopFer && nErrors >= 1, break; end
            end
        end

        function value = stat_field(s, names)
            value = 0;
            for ii = 1:numel(names)
                if isfield(s, names{ii})
                    value = s.(names{ii}); return;
                end
            end
        end

        function fer = sphere_packing_bound_fer(n, k, ebn0dBList)
            fer = zeros(size(ebn0dBList));
            dMin = max(3, min(2*floor((n-k)/8)+1, n-k+1));
            R = k/n;
            for ii = 1:numel(ebn0dBList)
                gamma = 10^(ebn0dBList(ii)/10);
                pb = 0.5*erfc(sqrt(2*R*dMin*gamma)/sqrt(2));
                fer(ii) = min(1, k*pb);
            end
        end
    end
end
