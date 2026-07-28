function coded_bits = cascade_encode(codec, msg_symbols)
    % RS outer-encode -> LSB-first bit serialize -> zero-pad -> BCH inner-encode each block.
    % Port of hc_src/cascade_v3.py CascadeV3Codec.encode.
    cfg = codec.cfg;
    assert(numel(msg_symbols) == cfg.k_rs, 'cascade_encode: bad message length');
    c_rs = shortened_rs_encode(codec.rs, msg_symbols);   % n_rs symbols

    rs_bits = zeros(1, codec.rs_bits);
    for i = 0:cfg.n_rs-1
        si = c_rs(i+1);
        for b = 0:cfg.m_rs-1
            rs_bits(i*cfg.m_rs + b + 1) = bitand(bitshift(si, -b), 1);
        end
    end
    if codec.n_pad_bits > 0
        rs_bits = [rs_bits, zeros(1, codec.n_pad_bits)];
    end

    coded_bits = zeros(1, codec.n_coded_bits);
    for b = 0:codec.n_bch_blocks-1
        block = rs_bits(b*codec.bch.k + 1 : (b+1)*codec.bch.k);
        coded_bits(b*codec.bch.n + 1 : (b+1)*codec.bch.n) = bch_encode(codec.bch, block);
    end
end
