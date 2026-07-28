function srs = shortened_rs_init(m, n_short, k_short)
    % Shortened narrow-sense RS(n_short,k_short) over GF(2^m).
    % Built from full RS(N=2^m-1, K_full) with the same parity count.
    % Port of hc_src/shortened_codes.py ShortenedRSCode.
    srs.m = m;
    srs.n = n_short;
    srs.k = k_short;
    srs.n_parity = n_short - k_short;
    assert(mod(srs.n_parity, 2) == 0, 'shortened_rs_init: RS parity must be even (=2t)');
    srs.t = srs.n_parity / 2;

    srs.N_full = bitshift(1, m) - 1;
    srs.K_full = srs.N_full - srs.n_parity;
    srs.shorten = srs.K_full - srs.k;
    assert(srs.shorten >= 0, 'shortened_rs_init: cannot shorten, k_short > K_full');
    assert(srs.N_full - srs.shorten == srs.n, 'shortened_rs_init: bad shortening');

    srs.full = rs_init(m, srs.K_full);
    srs.gf = srs.full.gf;
end
