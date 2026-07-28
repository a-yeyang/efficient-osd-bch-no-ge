function [msg_hat, ok_rs, n_bch_ok] = cascade_decode(codec, hard_bits)
    % BCH inner-decode each block -> reassemble RS symbols -> RS outer-decode.
    % Port of hc_src/cascade_v3.py CascadeV3Codec.decode.
    cfg = codec.cfg;

    recovered_bits = zeros(1, codec.n_bch_blocks * codec.bch.k);
    n_bch_ok = 0;
    for b = 0:codec.n_bch_blocks-1
        block = hard_bits(b*codec.bch.n + 1 : (b+1)*codec.bch.n);
        if strcmp(codec.bch_decoder, 'conv')
            [cw_hat, ok] = bch_decode_conventional_shortened(codec.bch, block);
        else
            [cw_hat, ok] = bch_decode_direct_shortened(codec.bch, block);
        end
        if ok
            n_bch_ok = n_bch_ok + 1;
        end
        recovered_bits(b*codec.bch.k + 1 : (b+1)*codec.bch.k) = bch_extract_message(codec.bch, cw_hat);
    end

    rs_bit_stream = recovered_bits(1:codec.rs_bits);
    r_rs = zeros(1, cfg.n_rs);
    for i = 0:cfg.n_rs-1
        s = 0;
        for b = 0:cfg.m_rs-1
            bitval = bitand(rs_bit_stream(i*cfg.m_rs + b + 1), 1);
            s = bitor(s, bitshift(bitval, b));
        end
        r_rs(i+1) = s;
    end

    [c_dec, ok_rs] = shortened_rs_decode(codec.rs, r_rs);
    if ok_rs
        msg_hat = shortened_rs_extract_message(codec.rs, c_dec);
    else
        msg_hat = r_rs(cfg.n_rs - cfg.k_rs + 1:end);
    end
end
