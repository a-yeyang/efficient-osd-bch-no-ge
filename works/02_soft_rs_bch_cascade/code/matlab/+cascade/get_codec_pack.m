function pack = get_codec_pack(cfg, method)
%GET_CODEC_PACK Build MATLAB function handles matching the Python simulator API.
method = lower(char(method));
switch method
    case {'pure_rs_bm', 'pure_rs_lccbr'}
        codec = cascade.PureRSCodec(cfg);
        if strcmp(method, 'pure_rs_bm')
            rsMethod = 'hard';
        else
            rsMethod = 'soft';
        end
        pack.encoder = @(message) codec.encode(message);
        pack.decoder = @(llr, counters) codec.decode(llr, rsMethod, counters);
        pack.effective_rate = codec.effective_rate;
        pack.codec = codec;
    case {'scheme_a', 'scheme_b'}
        codec = cascade.CascadedCodec(cfg);
        pack.encoder = @(message) codec.encode(message);
        pack.decoder = @(llr, counters) codec.decode(llr, method, counters);
        pack.effective_rate = codec.effective_rate;
        pack.codec = codec;
    otherwise
        error('cascade:Sim:UnknownMethod', 'Unknown method: %s', method);
end
end
