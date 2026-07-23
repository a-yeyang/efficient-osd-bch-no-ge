classdef RSCode < handle
    %RSCODE Narrow-sense primitive RS code with BM/Chien soft Chase decoding.
    properties (SetAccess = private)
        gf
        m
        n
        k
        d
        t
        alpha_pow
        g_poly
    end
    methods
        function obj = RSCode(m, k)
            obj.gf = cascade.GF(m);
            obj.m = double(m);
            obj.n = obj.gf.n;
            obj.k = double(k);
            assert(obj.k >= 1 && obj.k <= obj.n, 'cascade:RS:BadDimension', 'k must be in 1..n.');
            obj.d = obj.n - obj.k + 1;
            obj.t = floor((obj.d - 1) / 2);
            obj.alpha_pow = obj.gf.EXP(1:obj.n);
            obj.g_poly = 1;
            for ii = 1:(2 * obj.t)
                obj.g_poly = obj.gf.poly_mul(obj.g_poly, [obj.gf.alpha(ii), 1]);
            end
            assert(numel(obj.g_poly) - 1 == obj.n - obj.k, 'cascade:RS:GeneratorDegree', ...
                'RS generator degree did not match n-k.');
        end

        function codeword = encode_systematic(obj, message)
            message = double(message(:).');
            assert(numel(message) == obj.k, 'cascade:RS:MessageLength', 'Message has wrong length.');
            dividend = [zeros(1, obj.n - obj.k), message];
            [~, remainder] = obj.gf.poly_divmod(dividend, obj.g_poly);
            parity = [remainder, zeros(1, obj.n - obj.k - numel(remainder))];
            codeword = [parity, message];
        end

        function message = extract_message(obj, codeword)
            codeword = double(codeword(:).');
            message = codeword(obj.n - obj.k + 1:end);
        end

        function [decoded, ok] = bm_decode(obj, received, counters, cache)
            %BM_DECODE BM locator + Chien roots + Vandermonde magnitudes.
            %   Magnitudes are solved from the syndromes.  This is algebraically
            %   equivalent to Forney and avoids orientation ambiguity in a
            %   toolbox-free implementation while retaining the BM pipeline.
            if nargin < 3 || isempty(counters)
                counters = cascade.OpCounters();
            end
            if nargin < 4
                cache = [];
            end
            received = double(received(:).');
            assert(numel(received) == obj.n, 'cascade:RS:WordLength', 'Received word has wrong length.');
            nonzero = find(received ~= 0);
            if isempty(nonzero)
                decoded = received;
                ok = true;
                return;
            end
            syndromes = zeros(1, 2 * obj.t);
            for ii = 1:(2 * obj.t)
                value = 0;
                for jj = 1:numel(nonzero)
                    position0 = nonzero(jj) - 1;
                    exponentPosition = mod(ii * position0, obj.n) + 1;
                    if isempty(cache)
                        locator = obj.gf.alpha(exponentPosition - 1);
                    else
                        locator = cache.alpha_at(exponentPosition);
                    end
                    value = cascade.gfxor(value, obj.gf.mul(received(nonzero(jj)), locator));
                end
                syndromes(ii) = value;
            end
            counters.f2m = counters.f2m + numel(nonzero) * 2 * obj.t;
            if ~any(syndromes)
                decoded = received;
                ok = true;
                return;
            end

            [lambda, L] = cascade.berlekamp_massey(obj.gf, syndromes);
            nonzeroDegrees = find(lambda ~= 0) - 1;
            counters.f2m = counters.f2m + obj.n * numel(nonzeroDegrees);
            positions = cascade.chien_search(obj.gf, lambda, obj.n, cache);
            if numel(positions) ~= L || L > obj.t
                decoded = received;
                ok = false;
                return;
            end

            vandermonde = zeros(L, L);
            for row = 1:L
                for col = 1:L
                    exponentPosition = mod(row * (positions(col) - 1), obj.n) + 1;
                    if isempty(cache)
                        vandermonde(row, col) = obj.gf.alpha(exponentPosition - 1);
                    else
                        vandermonde(row, col) = cache.alpha_at(exponentPosition);
                    end
                end
            end
            [magnitudes, solveOK] = cascade.gf_solve(obj.gf, vandermonde, syndromes(1:L).');
            if ~solveOK
                decoded = received;
                ok = false;
                return;
            end
            decoded = received;
            for ii = 1:L
                decoded(positions(ii)) = cascade.gfxor(decoded(positions(ii)), magnitudes(ii));
            end
            % Reject a candidate if numerical/orientation errors left syndrome.
            valid = true;
            for ii = 1:(2 * obj.t)
                value = 0;
                nz = find(decoded ~= 0);
                for jj = 1:numel(nz)
                    exponentPosition = mod(ii * (nz(jj) - 1), obj.n) + 1;
                    if isempty(cache)
                        locator = obj.gf.alpha(exponentPosition - 1);
                    else
                        locator = cache.alpha_at(exponentPosition);
                    end
                    value = cascade.gfxor(value, obj.gf.mul(decoded(nz(jj)), locator));
                end
                if value ~= 0
                    valid = false;
                    break;
                end
            end
            if ~valid
                decoded = received;
                ok = false;
            else
                ok = true;
            end
        end

        function [decoded, ok] = lcc_br_decode(obj, received, reliability, eta, counters, cache)
            %LCC_BR_DECODE Simplified Chase-BM LCC-BR reference path.
            if nargin < 5 || isempty(counters)
                counters = cascade.OpCounters();
            end
            if nargin < 6
                cache = [];
            end
            received = double(received(:).');
            reliability = double(reliability(:).');
            eta = min(double(eta), obj.n);
            [~, order] = sort(reliability, 'ascend');
            leastReliable = order(1:eta);
            bestScore = -inf;
            bestCodeword = [];
            nTestVectors = 2^eta;
            for mask = 0:nTestVectors-1
                testWord = received;
                for bit = 1:eta
                    if bitget(uint32(mask), bit)
                        position = leastReliable(bit);
                        testWord(position) = cascade.gfxor(testWord(position), 1);
                    end
                end
                [candidate, candidateOK] = obj.bm_decode(testWord, counters, cache);
                if ~candidateOK
                    continue;
                end
                score = sum(reliability(candidate == received));
                if score > bestScore
                    bestScore = score;
                    bestCodeword = candidate;
                end
            end
            if isempty(bestCodeword)
                [candidate, candidateOK] = obj.bm_decode(received, counters, cache);
                if candidateOK
                    decoded = candidate;
                else
                    decoded = received;
                end
                ok = candidateOK;
                return;
            end
            counters.n_tvs = counters.n_tvs + nTestVectors;
            decoded = bestCodeword;
            ok = true;
        end
    end
end
