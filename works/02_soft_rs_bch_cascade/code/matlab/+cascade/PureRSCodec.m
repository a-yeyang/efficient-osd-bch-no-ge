classdef PureRSCodec < handle
    %PURERSCODEC Direct RS baseline used by the original comparison suite.
    properties (SetAccess = private)
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
            obj.rs = cascade.RSCode(cfg.m, cfg.k_rs);
            obj.rs_bits = cfg.m * cfg.n_rs;
            obj.n_coded_bits = obj.rs_bits;
            obj.n_info_bits = cfg.k_rs * cfg.m;
            obj.effective_rate = obj.n_info_bits / obj.n_coded_bits;
        end
        function bits = encode(obj, message_symbols)
            bits = cascade.symbols_to_bits(obj.rs.encode_systematic(message_symbols), obj.cfg.m);
        end
        function [message_hat, result] = decode(obj, llr, method, counters)
            if nargin < 3 || isempty(method)
                method = 'hard';
            end
            if nargin < 4 || isempty(counters)
                counters = cascade.OpCounters();
            end
            llr = double(llr(:).');
            hardBits = double(llr < 0);
            receivedSymbols = cascade.bits_to_symbols(hardBits, obj.cfg.m);
            switch lower(method)
                case 'hard'
                    [decodedRs, ok] = obj.rs.bm_decode(receivedSymbols, counters);
                case 'soft'
                    reliabilities = sum(reshape(abs(llr), obj.cfg.m, obj.cfg.n_rs), 1);
                    [decodedRs, ok] = obj.rs.lcc_br_decode(receivedSymbols, reliabilities, obj.cfg.lcc_eta, counters);
                otherwise
                    error('cascade:PureRS:UnknownMethod', 'Unknown PureRS decoder method: %s', method);
            end
            if ok
                message_hat = obj.rs.extract_message(decodedRs);
            else
                message_hat = receivedSymbols(obj.cfg.n_rs - obj.cfg.k_rs + 1:end);
            end
            result = struct('counters', counters, 'ok', logical(ok));
        end
    end
end
