classdef PureRSCodec < handle
%PURERSCODEC  RS-only baseline: no inner code (port of cascade_src/cascade.py
%   PureRSCodec). Coded bits = m*n_RS, rate = k_RS/n_RS.

    properties
        cfg
        rs
        rs_bits
        n_coded_bits
        n_info_bits
        effective_rate
    end

    methods
        function obj = PureRSCodec(cfg)
            obj.cfg = cfg;
            obj.rs = RSCode(cfg.m, cfg.k_rs);
            obj.rs_bits = cfg.m * cfg.n_rs();
            obj.n_coded_bits = obj.rs_bits;
            obj.n_info_bits = cfg.k_rs * cfg.m;
            obj.effective_rate = obj.n_info_bits / obj.n_coded_bits;
        end

        function bits = encode(obj, msg_symbols)
            cfg = obj.cfg;
            c_rs = obj.rs.encodeSystematic(msg_symbols);
            bits = zeros(1, cfg.m * cfg.n_rs());
            for i = 1:cfg.n_rs()
                s = c_rs(i);
                for b = 0:cfg.m-1
                    bits((i-1)*cfg.m + b + 1) = bitand(bitshift(s, -b), 1);
                end
            end
        end

        function [msg_hat, info] = decode(obj, llr, method, counters)
            % method: 'hard' (BM) or 'soft' (LCC-BR).
            if nargin < 4 || isempty(counters), counters = OpCounters(); end
            cfg = obj.cfg;
            hard_bits = double(llr < 0);
            r_rs = zeros(1, cfg.n_rs());
            for i = 0:cfg.n_rs()-1
                s = 0;
                for b = 0:cfg.m-1
                    s = bitor(s, bitshift(hard_bits(i*cfg.m + b + 1), b));
                end
                r_rs(i+1) = s;
            end

            if strcmp(method, 'hard')
                [c_dec, ok] = obj.rs.bmDecode(r_rs, counters);
            else  % soft
                abs_llr = abs(llr);
                reliability = sum(reshape(abs_llr, cfg.m, cfg.n_rs()), 1);
                [c_dec, ok] = obj.rs.lccBrDecode(r_rs, reliability, cfg.lcc_eta, counters);
            end

            if ok
                msg_hat = obj.rs.extractMessage(c_dec);
            else
                msg_hat = r_rs(cfg.n_rs() - cfg.k_rs + 1 : end);
            end
            info = struct('counters', counters, 'ok', ok);
        end
    end
end
