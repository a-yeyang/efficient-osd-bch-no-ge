classdef GF < handle
    %GF Finite-field arithmetic for GF(2^m), with integer polynomial symbols.
    %   The representation and primitive-polynomial table intentionally match
    %   Work 01 / the Python reference.  Elements are ordinary doubles in
    %   0..2^m-1 at the public API; addition is XOR.

    properties (SetAccess = private)
        m
        n
        prim
        EXP
        LOG
    end

    methods
        function obj = GF(m)
            primTable = [3, 7, 11, 19, 37, 67, 137, 285];
            validateattributes(m, {'numeric'}, {'scalar', 'integer', '>=', 1, '<=', 8});
            obj.m = double(m);
            obj.n = 2^obj.m - 1;
            obj.prim = uint32(primTable(obj.m));
            obj.EXP = zeros(1, 2 * obj.n + 2);
            obj.LOG = -ones(1, 2^obj.m);

            x = uint32(1);
            topBit = bitshift(uint32(1), obj.m);
            for ii = 0:obj.n-1
                obj.EXP(ii + 1) = double(x);
                obj.LOG(double(x) + 1) = ii;
                x = bitshift(x, 1);
                if bitand(x, topBit) ~= 0
                    x = bitxor(x, obj.prim);
                end
            end
            for ii = obj.n:(2 * obj.n + 1)
                obj.EXP(ii + 1) = obj.EXP(ii - obj.n + 1);
            end
        end

        function out = add(~, a, b)
            out = cascade.gfxor(a, b);
        end

        function out = mul(obj, a, b)
            if isscalar(a) && isscalar(b)
                if a == 0 || b == 0
                    out = 0;
                else
                    out = obj.EXP(obj.LOG(double(a) + 1) + obj.LOG(double(b) + 1) + 1);
                end
                return;
            end
            out = obj.vmul(a, b);
        end

        function out = inv(obj, a)
            if a == 0
                error('cascade:GF:ZeroInverse', 'Inverse of zero in GF(2^m).');
            end
            out = obj.EXP(obj.n - obj.LOG(double(a) + 1) + 1);
        end

        function out = div(obj, a, b)
            if b == 0
                error('cascade:GF:ZeroDivision', 'Division by zero in GF(2^m).');
            end
            if a == 0
                out = 0;
            else
                out = obj.EXP(mod(obj.LOG(double(a) + 1) - obj.LOG(double(b) + 1), obj.n) + 1);
            end
        end

        function out = pow(obj, a, exponent)
            if a == 0
                if exponent > 0
                    out = 0;
                else
                    out = 1;
                end
                return;
            end
            out = obj.EXP(mod(obj.LOG(double(a) + 1) * double(exponent), obj.n) + 1);
        end

        function out = vmul(obj, A, B)
            %VMUL Elementwise multiplication, with scalar expansion.
            if isscalar(A)
                A = repmat(A, size(B));
            elseif isscalar(B)
                B = repmat(B, size(A));
            elseif ~isequal(size(A), size(B))
                error('cascade:GF:SizeMismatch', 'GF operands must have compatible sizes.');
            end
            out = zeros(size(A));
            mask = A ~= 0 & B ~= 0;
            if any(mask, 'all')
                la = obj.LOG(double(A(mask)) + 1);
                lb = obj.LOG(double(B(mask)) + 1);
                out(mask) = obj.EXP(mod(la + lb, obj.n) + 1);
            end
        end

        function out = vinv(obj, A)
            if any(A == 0, 'all')
                error('cascade:GF:ZeroInverse', 'Inverse of zero in GF(2^m).');
            end
            out = obj.EXP(obj.n - obj.LOG(double(A) + 1) + 1);
        end

        function out = alpha(obj, exponent)
            %ALPHA Return alpha^exponent, accepting any integer exponent.
            out = obj.EXP(mod(double(exponent), obj.n) + 1);
        end

        function y = poly_eval(obj, coeffs, x)
            coeffs = double(coeffs(:).');
            y = 0;
            for ii = numel(coeffs):-1:1
                y = cascade.gfxor(obj.mul(y, x), coeffs(ii));
            end
        end

        function out = poly_mul(obj, a, b)
            a = double(a(:).');
            b = double(b(:).');
            out = zeros(1, numel(a) + numel(b) - 1);
            for ii = 1:numel(a)
                if a(ii) == 0
                    continue;
                end
                for jj = 1:numel(b)
                    if b(jj) ~= 0
                        out(ii + jj - 1) = cascade.gfxor(out(ii + jj - 1), obj.mul(a(ii), b(jj)));
                    end
                end
            end
        end

        function out = poly_add(~, a, b)
            a = double(a(:).');
            b = double(b(:).');
            out = zeros(1, max(numel(a), numel(b)));
            out(1:numel(a)) = a;
            out(1:numel(b)) = cascade.gfxor(out(1:numel(b)), b);
        end

        function [quot, rem] = poly_divmod(obj, numerator, denominator)
            numerator = double(numerator(:).');
            denominator = double(denominator(:).');
            while ~isempty(denominator) && denominator(end) == 0
                denominator(end) = [];
            end
            if isempty(denominator)
                error('cascade:GF:ZeroDivision', 'Zero polynomial divisor.');
            end
            qLen = max(0, numel(numerator) - numel(denominator) + 1);
            quot = zeros(1, qLen);
            rem = numerator;
            leadInv = obj.inv(denominator(end));
            for ii = qLen:-1:1
                idx = ii + numel(denominator) - 1;
                if idx > numel(rem)
                    continue;
                end
                coeff = obj.mul(rem(idx), leadInv);
                quot(ii) = coeff;
                if coeff ~= 0
                    for jj = 1:numel(denominator)
                        rem(ii + jj - 1) = cascade.gfxor(rem(ii + jj - 1), obj.mul(coeff, denominator(jj)));
                    end
                end
            end
            while ~isempty(rem) && rem(end) == 0
                rem(end) = [];
            end
        end
    end
end
