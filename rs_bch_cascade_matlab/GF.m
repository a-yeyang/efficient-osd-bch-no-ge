classdef GF
%GF  Finite field GF(2^m) built from a primitive polynomial.
%
%   Port of src/gf.py (GF class) from the Python reference implementation.
%
%   Elements are integers 0..2^m-1 (binary polynomial representation, bit i =
%   coefficient of alpha^i). EXP/LOG tables give O(1) multiply / inverse / pow;
%   addition is bitxor. Polynomials are little-endian row vectors: coeff of x^i
%   sits at index i+1 (MATLAB 1-based).
%
%   NOTE on indexing: field-element *values* and log *exponents* stay 0-based,
%   exactly as in Python. Only the raw table lookups add +1 for MATLAB. The
%   primitive polynomials match MATLAB's primpoly(m): primpoly(7)=137=0b10001001,
%   primpoly(8)=285=0b100011101, so this port is numerically consistent with
%   both the Python code and the Communications Toolbox.

    properties
        m       % GF exponent
        n       % 2^m - 1
        prim    % primitive polynomial (integer, bit i = coeff of x^i)
        EXP     % length 2n+2; EXP(i+1) = alpha^i, duplicated for mod-free lookup
        LOG     % length 2^m; LOG(x+1) = discrete log of x (0-based); LOG(1)=-1
    end

    methods
        function obj = GF(m)
            PRIM = containers.Map('KeyType','double','ValueType','double');
            PRIM(1)  = 0b11;          % x + 1
            PRIM(2)  = 0b111;         % x^2 + x + 1
            PRIM(3)  = 0b1011;        % x^3 + x + 1
            PRIM(4)  = 0b10011;       % x^4 + x + 1
            PRIM(5)  = 0b100101;      % x^5 + x^2 + 1
            PRIM(6)  = 0b1000011;     % x^6 + x + 1
            PRIM(7)  = 0b10001001;    % x^7 + x^3 + 1        (= 137, primpoly(7))
            PRIM(8)  = 0b100011101;   % x^8 + x^4+x^3+x^2+1  (= 285, primpoly(8))
            PRIM(9)  = 0b1000010001;  % x^9 + x^4 + 1
            PRIM(10) = 0b10000001001; % x^10 + x^3 + 1
            if ~isKey(PRIM, m)
                error('GF:noPrim', 'No primitive polynomial for m=%d', m);
            end
            obj.m = m;
            obj.n = bitshift(1, m) - 1;
            obj.prim = PRIM(m);

            EXP = zeros(1, 2*obj.n + 2);
            LOG = -1 * ones(1, bitshift(1, m));
            x = 1;
            for i = 0:obj.n-1
                EXP(i+1) = x;
                LOG(x+1) = i;
                x = bitshift(x, 1);
                if bitand(x, bitshift(1, m)) ~= 0
                    x = bitxor(x, obj.prim);
                end
            end
            % Duplicate so EXP(i+n) == EXP(i): mod-free lookups.
            for i = obj.n : (2*obj.n + 1)
                EXP(i+1) = EXP(i - obj.n + 1);
            end
            obj.EXP = EXP;
            obj.LOG = LOG;
        end

        % --- scalar ops ---------------------------------------------------
        function y = mul(obj, a, b)
            if a == 0 || b == 0
                y = 0;
            else
                y = obj.EXP(obj.LOG(a+1) + obj.LOG(b+1) + 1);
            end
        end

        function y = inv(obj, a)
            if a == 0
                error('GF:divZero', 'inverse of zero in GF(2^m)');
            end
            y = obj.EXP(obj.n - obj.LOG(a+1) + 1);
        end

        function y = div(obj, a, b)
            if a == 0
                y = 0; return;
            end
            if b == 0
                error('GF:divZero', 'divide by zero in GF(2^m)');
            end
            y = obj.EXP(obj.LOG(a+1) - obj.LOG(b+1) + obj.n + 1);
        end

        function y = powr(obj, a, e)
            % a^e in GF(2^m). (Named powr to avoid clashing with builtin pow.)
            if a == 0
                if e > 0, y = 0; else, y = 1; end
                return;
            end
            y = obj.EXP(mod(obj.LOG(a+1) * e, obj.n) + 1);
        end

        % --- vectorized table lookups (values/exponents stay 0-based) -----
        function v = expv(obj, idx0)
            % EXP at 0-based exponent indices idx0 (may exceed n; caller ensures < 2n+2)
            v = obj.EXP(idx0 + 1);
        end

        function l = logv(obj, vals)
            % discrete logs (0-based) of field-element values (must be nonzero)
            l = obj.LOG(vals + 1);
        end

        function C = vmul(obj, A, B)
            % element-wise GF multiply of equal-shape (or broadcastable) arrays
            A = double(A); B = double(B);
            C = zeros(size(A .* ones(size(B))));
            if isscalar(A), A = A * ones(size(C)); end
            if isscalar(B), B = B * ones(size(C)); end
            mask = (A ~= 0) & (B ~= 0);
            if any(mask(:))
                la = obj.LOG(A(mask) + 1);
                lb = obj.LOG(B(mask) + 1);
                C(mask) = obj.EXP(la + lb + 1);
            end
        end

        % --- polynomial ops over GF(2^m) (little-endian row vectors) ------
        function y = polyEval(obj, coeffs, x)
            % Horner: p(x) with coeffs(i+1) = coeff of x^i
            y = 0;
            for idx = numel(coeffs):-1:1
                y = bitxor(obj.mul(y, x), coeffs(idx));
            end
        end

        function out = polyMul(obj, a, b)
            a = a(:).'; b = b(:).';
            out = zeros(1, numel(a) + numel(b) - 1);
            for i = 1:numel(a)
                if a(i) == 0, continue; end
                for j = 1:numel(b)
                    if b(j) == 0, continue; end
                    out(i+j-1) = bitxor(out(i+j-1), obj.mul(a(i), b(j)));
                end
            end
        end

        function out = polyAdd(~, a, b)
            a = a(:).'; b = b(:).';
            L = max(numel(a), numel(b));
            a(end+1:L) = 0; b(end+1:L) = 0;
            out = bitxor(a, b);
        end

        function [quot, rem] = polyDivmod(obj, num, den)
            num = num(:).'; den = den(:).';
            while ~isempty(den) && den(end) == 0, den(end) = []; end
            if isempty(den), error('GF:divZero', 'zero divisor'); end
            lead_inv = obj.inv(den(end));
            ld = numel(den);
            qlen = max(0, numel(num) - ld + 1);
            quot = zeros(1, qlen);
            rem = num;
            for i = qlen-1:-1:0    % i is 0-based (matches Python)
                if numel(rem) <= i + ld - 1
                    continue;
                end
                coeff = obj.mul(rem(i + ld - 1 + 1), lead_inv);
                quot(i+1) = coeff;
                for j = 0:ld-1
                    rem(i + j + 1) = bitxor(rem(i + j + 1), obj.mul(coeff, den(j+1)));
                end
            end
            while ~isempty(rem) && rem(end) == 0, rem(end) = []; end
        end

        % --- BCH generator polynomial (binary, deg = n-k) -----------------
        function g = bchGeneratorPoly(obj, t)
            % Product of minimal polynomials of alpha^1,alpha^3,...,alpha^{2t-1}.
            seen = [];
            g = 1;   % constant polynomial 1
            for i = 1:(2*t)
                % cyclotomic coset of i mod n
                coset = [];
                j = i;
                while ~ismember(j, coset)
                    coset(end+1) = j; %#ok<AGROW>
                    j = mod(j * 2, obj.n);
                end
                rep = min(coset);
                if ismember(rep, seen), continue; end
                seen(end+1) = rep; %#ok<AGROW>
                % minimal polynomial = prod_{s in coset} (x - alpha^s), char 2
                mpoly = 1;
                for s = coset
                    root = obj.EXP(s + 1);
                    mpoly = obj.polyMul(mpoly, [root, 1]);
                end
                % coefficients are guaranteed binary
                mpoly = mod(mpoly, 2);
                g = obj.polyMul(g, mpoly);
                g = mod(g, 2);
            end
        end
    end
end
