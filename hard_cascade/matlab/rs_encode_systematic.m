function c = rs_encode_systematic(rs, msg)
    % Systematic RS encode: c = [parity(n-k), msg(k)].
    assert(numel(msg) == rs.k, 'rs_encode_systematic: bad message length');
    gf = rs.gf;
    nk = rs.n - rs.k;
    dividend = [zeros(1, nk), double(msg(:).')];
    [~, rem] = gf_polydivmod(gf, dividend, rs.g_poly);
    parity = zeros(1, nk);
    parity(1:numel(rem)) = rem;
    c = zeros(1, rs.n);
    c(1:nk) = parity;
    c(nk+1:end) = msg(:).';
end
