function rs = rs_init(m, k)
    % Narrow-sense primitive RS(n,k) over GF(2^m), n = 2^m - 1.
    % Port of cascade_src/rs_code.py RSCode.__init__ / rs_bch_cascade_matlab/RSCode.m
    rs.gf = GF_init(m);
    rs.m = m;
    rs.n = rs.gf.n;
    rs.k = k;
    rs.d = rs.n - rs.k + 1;
    rs.t = floor((rs.d - 1) / 2);

    % g(x) = prod_{i=1..2t} (x - alpha^i); char 2 => minus = plus.
    g = 1;
    for i = 1:(2*rs.t)
        root = rs.gf.EXP(i + 1);   % alpha^i (0-based exponent i)
        g = gf_polymul(rs.gf, g, [root, 1]);
    end
    rs.g_poly = g;
    assert(numel(rs.g_poly) - 1 == rs.n - rs.k, 'rs_init: g_poly degree mismatch');
end
