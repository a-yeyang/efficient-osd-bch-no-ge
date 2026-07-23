classdef CascadedCodec < handle
    %CASCADEDCODEC RS outer code + BCH/LLOSD inner code and both schemes.
    properties (SetAccess = private)
        cfg
        rs
        bch
        rs_bits
        n_pad_bits
        n_bch_blocks
        n_coded_bits
        n_info_bits
        effective_rate
        bch_pivots
        bch_recover_matrix
    end
    methods
        function obj = CascadedCodec(cfg)
            obj.cfg = cfg;
            obj.rs = cascade.RSCode(cfg.m, cfg.k_rs);
            obj.bch = cascade.BCHCode(cfg.m, cfg.t_bch);
            obj.rs_bits = cfg.m * cfg.n_rs;
            remainder = mod(obj.rs_bits, obj.bch.k);
            if remainder == 0
                obj.n_pad_bits = 0;
            else
                obj.n_pad_bits = obj.bch.k - remainder;
            end
            obj.n_bch_blocks = (obj.rs_bits + obj.n_pad_bits) / obj.bch.k;
            obj.n_coded_bits = obj.n_bch_blocks * obj.bch.n;
            obj.n_info_bits = cfg.k_rs * cfg.m;
            obj.effective_rate = obj.n_info_bits / obj.n_coded_bits;
            [obj.bch_pivots, obj.bch_recover_matrix] = cascade.gf2_independent_columns(obj.bch.G);
        end

        function bits = encode(obj, message_symbols)
            message_symbols = double(message_symbols(:).');
            assert(numel(message_symbols) == obj.cfg.k_rs, 'cascade:Cascade:MessageLength', 'Message has wrong length.');
            rsCodeword = obj.rs.encode_systematic(message_symbols);
            rsBits = cascade.symbols_to_bits(rsCodeword, obj.cfg.m);
            if obj.n_pad_bits > 0
                rsBits = [rsBits, zeros(1, obj.n_pad_bits)];
            end
            bits = zeros(1, obj.n_coded_bits);
            for block = 1:obj.n_bch_blocks
                source = rsBits((block - 1) * obj.bch.k + 1:block * obj.bch.k);
                bits((block - 1) * obj.bch.n + 1:block * obj.bch.n) = obj.bch.encode(source);
            end
        end

        function [message_hat, result] = decode(obj, llr, method, counters)
            if nargin < 3 || isempty(method)
                method = 'scheme_a';
            end
            if nargin < 4 || isempty(counters)
                counters = cascade.OpCounters();
            end
            llr = double(llr(:).');
            assert(numel(llr) == obj.n_coded_bits, 'cascade:Cascade:WordLength', 'LLR vector has wrong length.');
            switch lower(method)
                case 'scheme_a'
                    [message_hat, result] = obj.decode_scheme_a(llr, counters);
                case 'scheme_b'
                    [message_hat, result] = cascade.cascade_scheme_b_decode(obj, llr, [], counters);
                otherwise
                    error('cascade:Cascade:UnknownMethod', 'Unknown cascade decoder method: %s', method);
            end
        end

        function [message_hat, result] = decode_scheme_a(obj, llr, counters, cache)
            %DECODE_SCHEME_A LLOSD inner + Chase/BM LCC-BR outer, no sharing.
            if nargin < 4
                cache = [];
            end
            decodedBits = zeros(1, obj.n_bch_blocks * obj.bch.k);
            blockReliability = zeros(1, obj.n_bch_blocks);
            for block = 1:obj.n_bch_blocks
                blockLlr = llr((block - 1) * obj.bch.n + 1:block * obj.bch.n);
                [decodedCodeword, innerStats] = cascade.llosd_fast(obj.bch, blockLlr, obj.cfg.llosd_tau, true, true, cache);
                % Unlike the earlier prototype, account for every inner
                % LLOSD invocation.  This makes Monte-Carlo/KPI numbers cover
                % the full cascade rather than only outer RS Chase-BM work.
                counters.add(innerStats.counters);
                decodedBits((block - 1) * obj.bch.k + 1:block * obj.bch.k) = obj.bch_msg_recover(decodedCodeword);
                blockReliability(block) = sum(abs(blockLlr));
            end
            receivedSymbols = obj.bit_stream_to_rs_symbols(decodedBits);
            symbolReliability = zeros(1, obj.cfg.n_rs);
            for symbol = 1:obj.cfg.n_rs
                bitIndices = (symbol - 1) * obj.cfg.m + 1:symbol * obj.cfg.m;
                blockIds = unique(floor((bitIndices - 1) / obj.bch.k) + 1);
                symbolReliability(symbol) = mean(blockReliability(blockIds));
            end
            [decodedRs, ok] = obj.rs.lcc_br_decode(receivedSymbols, symbolReliability, obj.cfg.lcc_eta, counters, cache);
            if ok
                message_hat = obj.rs.extract_message(decodedRs);
            else
                message_hat = receivedSymbols(obj.cfg.n_rs - obj.cfg.k_rs + 1:end);
            end
            result = struct('counters', counters, 'ok', logical(ok));
        end

        function symbols = bit_stream_to_rs_symbols(obj, bits)
            bits = double(bits(:).');
            useful = bits(1:obj.rs_bits);
            symbols = cascade.bits_to_symbols(useful, obj.cfg.m);
        end

        function message = bch_msg_recover(obj, codeword)
            codeword = double(codeword(:).');
            message = mod(codeword(obj.bch_pivots) * obj.bch_recover_matrix, 2);
        end
    end
end
