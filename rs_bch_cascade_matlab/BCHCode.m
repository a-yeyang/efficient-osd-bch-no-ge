classdef BCHCode
%BCHCODE  Primitive narrow-sense binary BCH code of length n = 2^m - 1.
%
%   Port of src/bch.py (BCHCode class). NON-systematic generator matrix
%   G (k x n): row i = g(x) * x^i, taken mod 2, exactly as in the reference.
%   Parity-check H is built from alpha^{ij} expanded to binary and row-reduced.
%   Codewords are 0/1 row vectors. bmDecode gives hard-decision decoding.

    properties
        m
        t
        gf
        n
        g_poly      % binary generator polynomial (little-endian, coeffs 0/1)
        k
        d_design    % 2t + 1
        G           % k x n binary generator matrix
        H           % (n-k) x n binary parity-check matrix
        % msg-recover cache (for cascade): pivots + inverse of G(:,pivots)
        pivots
        msgRecoverMatrix
    end

    methods
        function obj = BCHCode(m, t)
            obj.m = m;
            obj.t = t;
            obj.gf = GF(m);
            obj.n = obj.gf.n;
            obj.g_poly = obj.gf.bchGeneratorPoly(t);   % coeffs in {0,1}
            deg = numel(obj.g_poly) - 1;
            while deg > 0 && obj.g_poly(deg + 1) == 0, deg = deg - 1; end
            obj.k = obj.n - deg;
            obj.d_design = 2*t + 1;

            % Generator matrix: row i = g(x) * x^i.
            g = double(obj.g_poly(:).');
            assert(numel(g) == obj.n - obj.k + 1);
            G = zeros(obj.k, obj.n);
            for i = 0:obj.k-1
                G(i+1, (i+1):(i+numel(g))) = g;
            end
            obj.G = G;

            % Parity-check H: alpha^{ij} for i=1..2t expanded to m binary rows.
            gf = obj.gf;
            H_ext = zeros(obj.m * (2*t), obj.n);
            for j = 0:obj.n-1
                for i = 1:(2*t)
                    e = gf.EXP(mod(i*j, gf.n) + 1);   % alpha^{ij}
                    for b = 0:obj.m-1
                        H_ext((i-1)*obj.m + b + 1, j + 1) = bitand(bitshift(e, -b), 1);
                    end
                end
            end
            obj.H = BCHCode.rowReduceBinary(H_ext);
            assert(all(all(mod(obj.G * obj.H.', 2) == 0)), 'G H^T != 0');
        end

        function c = encode(obj, msg)
            c = mod(double(msg(:).') * obj.G, 2);
        end

        % -----------------------------------------------------------------
        function [c, ok] = bmDecode(obj, r_hard)
            % Berlekamp-Massey hard-decision decode for binary BCH.
            gf = obj.gf; t = obj.t;
            r_hard = double(r_hard(:).');
            S = zeros(1, 2*t + 1);
            idx = find(r_hard ~= 0) - 1;   % 0-based error-candidate positions
            for i = 1:(2*t)
                s = 0;
                for jj = 1:numel(idx)
                    s = bitxor(s, gf.EXP(mod(i*idx(jj), gf.n) + 1));
                end
                S(i+1) = s;
            end
            if ~any(S(2:end))
                c = r_hard; ok = true; return;
            end

            L = 0; Lam = 1; B = 1; b = 1; m_shift = 1;
            for nn = 1:(2*t)
                delta = S(nn+1);
                for i = 1:L
                    if (i+1) <= numel(Lam) && Lam(i+1) ~= 0
                        delta = bitxor(delta, gf.mul(Lam(i+1), S(nn - i + 1)));
                    end
                end
                if delta == 0
                    m_shift = m_shift + 1;
                else
                    coef = gf.div(delta, b);
                    xmB = [zeros(1, m_shift), B];
                    new_len = max(numel(Lam), numel(xmB));
                    T = [Lam, zeros(1, new_len - numel(Lam))];
                    xmB = [xmB, zeros(1, new_len - numel(xmB))];
                    for i = 1:new_len
                        T(i) = bitxor(T(i), gf.mul(coef, xmB(i)));
                    end
                    if 2*L <= nn - 1
                        L_new = nn - L; B = Lam; b = delta; Lam = T; L = L_new; m_shift = 1;
                    else
                        Lam = T; m_shift = m_shift + 1;
                    end
                end
            end

            err_positions = [];
            for i = 0:obj.n-1
                val = 0;
                for jj = 0:numel(Lam)-1
                    if Lam(jj+1) == 0, continue; end
                    val = bitxor(val, gf.mul(Lam(jj+1), gf.EXP(mod((obj.n - i)*jj, gf.n) + 1)));
                end
                if val == 0, err_positions(end+1) = i; end %#ok<AGROW>
            end
            if numel(err_positions) ~= L || L > t
                c = r_hard; ok = false; return;
            end
            c = r_hard; ok = true;
            for pidx = 1:numel(err_positions)
                p = err_positions(pidx);
                c(p+1) = bitxor(c(p+1), 1);
            end
        end

        % -----------------------------------------------------------------
        function obj = buildMsgRecover(obj)
            % Precompute pivots + inverse to recover k-bit message from a
            % non-systematic BCH codeword c = msg @ G.
            G = mod(obj.G, 2);
            k = obj.k; n = obj.n;
            pivots = [];
            for col = 0:n-1
                if numel(pivots) == k, break; end
                piv_row = [];
                for r = numel(pivots):k-1
                    if G(r+1, col+1) ~= 0, piv_row = r; break; end
                end
                if isempty(piv_row), continue; end
                pivots(end+1) = col; %#ok<AGROW>
                lead = numel(pivots) - 1;
                if piv_row ~= lead
                    G([lead+1, piv_row+1], :) = G([piv_row+1, lead+1], :);
                end
                for r = 0:k-1
                    if r ~= lead && G(r+1, col+1) ~= 0
                        G(r+1, :) = bitxor(G(r+1, :), G(lead+1, :));
                    end
                end
            end
            assert(numel(pivots) == k, 'BCH G rank < k');

            % inverse of G_orig(:, pivots) mod 2 via Gauss-Jordan
            G_orig = mod(obj.G, 2);
            Gp = G_orig(:, pivots + 1);
            inv = eye(k);
            for i = 0:k-1
                piv = [];
                for r = i:k-1
                    if Gp(r+1, i+1) ~= 0, piv = r; break; end
                end
                assert(~isempty(piv));
                if piv ~= i
                    Gp([i+1, piv+1], :) = Gp([piv+1, i+1], :);
                    inv([i+1, piv+1], :) = inv([piv+1, i+1], :);
                end
                for r = 0:k-1
                    if r ~= i && Gp(r+1, i+1) ~= 0
                        Gp(r+1, :) = bitxor(Gp(r+1, :), Gp(i+1, :));
                        inv(r+1, :) = bitxor(inv(r+1, :), inv(i+1, :));
                    end
                end
            end
            obj.pivots = pivots;
            obj.msgRecoverMatrix = inv;
        end

        function msg = recoverMsg(obj, c_hat)
            % Recover k-bit message; requires buildMsgRecover to have run.
            c_hat = double(c_hat(:).');
            msg = mod(c_hat(obj.pivots + 1) * obj.msgRecoverMatrix, 2);
        end
    end

    methods (Static)
        function A = rowReduceBinary(M)
            A = mod(M, 2);
            [rows, cols] = size(A);
            r = 0;
            for c = 0:cols-1
                if r >= rows, break; end
                piv = [];
                for i = r:rows-1
                    if A(i+1, c+1) ~= 0, piv = i; break; end
                end
                if isempty(piv), continue; end
                if piv ~= r
                    A([r+1, piv+1], :) = A([piv+1, r+1], :);
                end
                for i = 0:rows-1
                    if i ~= r && A(i+1, c+1) ~= 0
                        A(i+1, :) = bitxor(A(i+1, :), A(r+1, :));
                    end
                end
                r = r + 1;
            end
            keep = any(A, 2);
            A = A(keep, :);
        end
    end
end
