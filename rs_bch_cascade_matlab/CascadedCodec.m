classdef CascadedCodec < handle
%CASCADEDCODEC  RS+BCH cascade encoder/decoder (port of cascade_src/cascade.py
%   CascadedCodec). Outer RS (symbol-level), inner BCH (bit-level).
%
%   Encoding chain: msg symbols -> RS systematic encode -> bits (LSB first,
%   m bits/symbol) -> zero-pad to a multiple of BCH.k -> BCH-encode each block
%   -> concatenated coded bit stream.
%
%   Decoding chain (decode): per-bit LLR -> per BCH block: inner soft decode
%   (LLOSD or OSD, selectable) -> recover k-bit BCH message -> reassemble RS
%   bit stream (strip padding) -> RS symbols -> RS LCC-BR outer decode.

    properties
        cfg
        rs
        bch
        rs_bits
        n_pad_bits
        n_bch_blocks
        n_coded_bits
        n_info_bits
        effective_rate
    end

    methods
        function obj = CascadedCodec(cfg)
            obj.cfg = cfg;
            obj.rs = RSCode(cfg.m, cfg.k_rs);
            obj.bch = BCHCode(cfg.m, cfg.t_bch);
            obj.bch = obj.bch.buildMsgRecover();

            obj.rs_bits = cfg.m * cfg.n_rs();
            if mod(obj.rs_bits, obj.bch.k) ~= 0
                obj.n_pad_bits = obj.bch.k - mod(obj.rs_bits, obj.bch.k);
            else
                obj.n_pad_bits = 0;
            end
            obj.n_bch_blocks = (obj.rs_bits + obj.n_pad_bits) / obj.bch.k;
            obj.n_coded_bits = obj.n_bch_blocks * obj.bch.n;
            obj.n_info_bits = cfg.k_rs * cfg.m;
            obj.effective_rate = obj.n_info_bits / obj.n_coded_bits;
        end

        % -----------------------------------------------------------------
        function bits = encode(obj, msg_symbols)
            cfg = obj.cfg;
            assert(numel(msg_symbols) == cfg.k_rs);

            c_rs = obj.rs.encodeSystematic(msg_symbols);   % n_rs symbols

            % symbols -> bits, m bits per symbol, LSB first
            rs_bits_v = zeros(1, cfg.m * cfg.n_rs());
            for i = 1:cfg.n_rs()
                s = c_rs(i);
                for b = 0:cfg.m-1
                    rs_bits_v((i-1)*cfg.m + b + 1) = bitand(bitshift(s, -b), 1);
                end
            end

            if obj.n_pad_bits > 0
                rs_bits_padded = [rs_bits_v, zeros(1, obj.n_pad_bits)];
            else
                rs_bits_padded = rs_bits_v;
            end

            bits = zeros(1, obj.n_coded_bits);
            for b = 0:obj.n_bch_blocks-1
                block = rs_bits_padded(b*obj.bch.k + 1 : (b+1)*obj.bch.k);
                cw = obj.bch.encode(block);
                bits(b*obj.bch.n + 1 : (b+1)*obj.bch.n) = cw;
            end
        end

        % -----------------------------------------------------------------
        function [msg_hat, info] = decode(obj, llr, innerMethod, counters)
            % innerMethod: 'llosd' (Lagrange, the paper's method) or
            %              'osd' (traditional Gaussian-elimination contrast).
            if nargin < 4 || isempty(counters), counters = OpCounters(); end
            cfg = obj.cfg;
            bch = obj.bch;
            assert(numel(llr) == obj.n_coded_bits);

            bch_decoded_bits = zeros(1, obj.n_bch_blocks * bch.k);
            block_reliability = zeros(1, obj.n_bch_blocks);

            for b = 0:obj.n_bch_blocks-1
                block_llr = llr(b*bch.n + 1 : (b+1)*bch.n);
                switch innerMethod
                    case 'llosd'
                        [c_hat, ~] = llosd_decode(bch, block_llr, cfg.llosd_tau, counters);
                    case 'osd'
                        [c_hat, ~] = osd_decode(bch, block_llr, cfg.llosd_tau, counters);
                    otherwise
                        error('CascadedCodec:innerMethod', 'unknown inner method %s', innerMethod);
                end
                msg_bits = bch.recoverMsg(c_hat);
                bch_decoded_bits(b*bch.k + 1 : (b+1)*bch.k) = msg_bits;
                block_reliability(b+1) = sum(abs(block_llr));
            end

            r_rs = obj.bitStreamToRsSymbols(bch_decoded_bits);

            symbol_reliability = zeros(1, cfg.n_rs());
            for i = 0:cfg.n_rs()-1
                bit_range = (i*cfg.m) : (i*cfg.m + cfg.m - 1);   % 0-based bit indices
                block_ids = floor(bit_range / bch.k);            % 0-based block ids
                symbol_reliability(i+1) = mean(block_reliability(block_ids + 1));
            end

            [c_dec, ok] = obj.rs.lccBrDecode(r_rs, symbol_reliability, cfg.lcc_eta, counters);
            if ok
                msg_hat = obj.rs.extractMessage(c_dec);
            else
                msg_hat = r_rs(cfg.n_rs() - cfg.k_rs + 1 : end);
            end
            info = struct('counters', counters, 'ok', ok);
        end

        % -----------------------------------------------------------------
        function r_rs = bitStreamToRsSymbols(obj, bits)
            cfg = obj.cfg;
            useful_bits = bits(1:obj.rs_bits);
            r_rs = zeros(1, cfg.n_rs());
            for i = 0:cfg.n_rs()-1
                s = 0;
                for b = 0:cfg.m-1
                    bitv = mod(round(useful_bits(i*cfg.m + b + 1)), 2);
                    s = bitor(s, bitshift(bitv, b));
                end
                r_rs(i+1) = s;
            end
        end
    end
end
