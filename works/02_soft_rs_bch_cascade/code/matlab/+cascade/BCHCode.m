classdef BCHCode < handle
    %BCHCODE Primitive narrow-sense BCH code with a nonsystematic G matrix.
    properties (SetAccess = private)
        m
        t
        gf
        n
        g_poly
        k
        d_design
        G
        H
    end
    methods
        function obj = BCHCode(m, t)
            obj.m = double(m);
            obj.t = double(t);
            obj.gf = cascade.GF(m);
            obj.n = obj.gf.n;
            obj.g_poly = cascade.bch_generator_poly(obj.gf, obj.t);
            obj.k = cascade.bch_dimension(obj.g_poly, obj.n);
            obj.d_design = 2 * obj.t + 1;

            obj.G = zeros(obj.k, obj.n);
            for ii = 1:obj.k
                obj.G(ii, ii:ii + numel(obj.g_poly) - 1) = obj.g_poly;
            end

            hExtended = zeros(obj.m * 2 * obj.t, obj.n);
            for jj = 0:obj.n-1
                for ii = 1:(2 * obj.t)
                    element = obj.gf.alpha(mod(ii * jj, obj.n));
                    for bit = 0:obj.m-1
                        hExtended((ii - 1) * obj.m + bit + 1, jj + 1) = bitget(uint32(element), bit + 1);
                    end
                end
            end
            obj.H = cascade.row_reduce_binary(hExtended);
            if ~all(mod(obj.G * obj.H.', 2) == 0, 'all')
                error('cascade:BCH:ParityCheckFailure', 'Generated BCH G and H are incompatible.');
            end
        end

        function codeword = encode(obj, message)
            message = mod(double(message(:).'), 2);
            assert(numel(message) == obj.k, 'cascade:BCH:MessageLength', 'Message has wrong length.');
            codeword = mod(message * obj.G, 2);
        end

        function [decoded, ok] = bm_decode(obj, received)
            %BM_DECODE Binary bounded-distance BM/Chien decoder.
            received = mod(double(received(:).'), 2);
            assert(numel(received) == obj.n, 'cascade:BCH:WordLength', 'Received word has wrong length.');
            gf = obj.gf;
            syndromes = zeros(1, 2 * obj.t);
            nonzero = find(received ~= 0) - 1;
            for ii = 1:(2 * obj.t)
                value = 0;
                for jj = 1:numel(nonzero)
                    value = cascade.gfxor(value, gf.alpha(mod(ii * nonzero(jj), gf.n)));
                end
                syndromes(ii) = value;
            end
            if ~any(syndromes)
                decoded = received;
                ok = true;
                return;
            end

            [lambda, L] = cascade.berlekamp_massey(gf, syndromes);
            errorPositions = cascade.chien_search(gf, lambda, obj.n);
            if numel(errorPositions) ~= L || L > obj.t
                decoded = received;
                ok = false;
                return;
            end
            decoded = received;
            decoded(errorPositions) = mod(decoded(errorPositions) + 1, 2);
            ok = true;
        end
    end
end
