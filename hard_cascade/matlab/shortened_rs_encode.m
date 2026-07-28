function c_short = shortened_rs_encode(srs, msg)
    assert(numel(msg) == srs.k, 'shortened_rs_encode: bad message length');
    full_msg = [zeros(1, srs.shorten), double(msg(:).')];
    c_full = rs_encode_systematic(srs.full, full_msg);
    % c_full = [ parity(n_parity) | 0*shorten | real_msg ]
    parity = c_full(1:srs.n_parity);
    real = c_full(srs.n_parity + srs.shorten + 1:end);
    c_short = [parity, real];
end
