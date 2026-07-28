function cw_short = bch_encode(bch, msg)
    assert(numel(msg) == bch.k, 'bch_encode: bad message length');
    full_msg = [zeros(1, bch.shorten), double(msg(:).')];
    c_full = mod(full_msg * bch.G, 2);   % [ parity | 0*shorten | real_msg ]
    parity = c_full(1:bch.n_parity);
    real = c_full(bch.n_parity + bch.shorten + 1:end);
    cw_short = [parity, real];
end
