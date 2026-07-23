classdef LagrangeCache < handle
    %LAGRANGECACHE Shared GF/Lagrange algebra cache for Scheme B.
    %   The reference implementation uses it primarily for accounting; this
    %   class additionally exposes the actual denominator and basis cache so
    %   it can be used by future LLOSD/LCC-BR refinements without API changes.
    properties (SetAccess = private)
        gf
        n
        alpha_pow
        pairwise_diff
        denom_cache
        ops_saved = 0
        alpha_hits = 0
        basis_queries = 0
    end
    methods
        function obj = LagrangeCache(gf, n)
            obj.gf = gf;
            obj.n = double(n);
            obj.alpha_pow = gf.EXP(mod(0:obj.n-1, gf.n) + 1);
            obj.pairwise_diff = [];
            obj.denom_cache = containers.Map('KeyType', 'char', 'ValueType', 'double');
        end

        function table = build_pairwise_diff(obj)
            if isempty(obj.pairwise_diff)
                left = repmat(obj.alpha_pow(:), 1, obj.n);
                right = repmat(obj.alpha_pow, obj.n, 1);
                obj.pairwise_diff = cascade.gfxor(left, right);
            end
            table = obj.pairwise_diff;
        end

        function value = denominator_product(obj, position, support)
            support = double(support(:).');
            key = sprintf('%d|%s', position, sprintf('%d,', sort(support)));
            if isKey(obj.denom_cache, key)
                % The denominator was already formed for another parity-point
                % evaluation.  Count one canonical GF-operation unit per
                % reuse; this matches the operation-accounting granularity in
                % LLOSD rather than charging the MATLAB loop implementation.
                obj.ops_saved = obj.ops_saved + 1;
                value = obj.denom_cache(key);
                return;
            end
            value = 1;
            for ii = 1:numel(support)
                if support(ii) ~= position
                    difference = obj.difference_at(position, support(ii));
                    value = obj.gf.mul(value, difference);
                end
            end
            obj.denom_cache(key) = value;
        end

        function value = lagrange_basis(obj, basis_position, support, eval_position)
            support = double(support(:).');
            obj.basis_queries = obj.basis_queries + 1;
            numerator = 1;
            for ii = 1:numel(support)
                if support(ii) ~= basis_position
                    difference = obj.difference_at(eval_position, support(ii));
                    numerator = obj.gf.mul(numerator, difference);
                end
            end
            denominator = obj.denominator_product(basis_position, support);
            value = obj.gf.div(numerator, denominator);
        end

        function table = lagrange_basis_table(obj, info_positions, eval_positions)
            info_positions = double(info_positions(:).');
            eval_positions = double(eval_positions(:).');
            table = zeros(numel(eval_positions), numel(info_positions));
            for ii = 1:numel(eval_positions)
                for jj = 1:numel(info_positions)
                    table(ii, jj) = obj.lagrange_basis(info_positions(jj), info_positions, eval_positions(ii));
                end
            end
        end

        function value = alpha_at(obj, position)
            %ALPHA_AT Return a cached locator at a one-based code position.
            position = mod(double(position) - 1, obj.n) + 1;
            obj.alpha_hits = obj.alpha_hits + 1;
            value = obj.alpha_pow(position);
        end

        function out = stats(obj)
            out = struct('cache_size', obj.denom_cache.Count, 'ops_saved', obj.ops_saved, ...
                'alpha_hits', obj.alpha_hits, 'basis_queries', obj.basis_queries, ...
                'pairwise_table_built', ~isempty(obj.pairwise_diff));
        end
    end

    methods (Access = private)
        function value = difference_at(obj, leftPosition, rightPosition)
            %DIFFERENCE_AT alpha^(left-1) + alpha^(right-1), using the shared
            %pairwise table after it has been initialized by Scheme B.
            if isempty(obj.pairwise_diff)
                value = cascade.gfxor(obj.alpha_pow(leftPosition), obj.alpha_pow(rightPosition));
            else
                value = obj.pairwise_diff(leftPosition, rightPosition);
            end
        end
    end
end
