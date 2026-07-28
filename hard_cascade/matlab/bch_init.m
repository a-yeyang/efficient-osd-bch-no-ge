function bch = bch_init(m, t, n_short, k_short)
    % Binary primitive narrow-sense BCH code, t in {1,2}, optional shortening.
    % Port of hc_src/shortened_codes.py BinaryBCH.__init__.
    assert(t == 1 || t == 2, 'bch_init: only t=1 and t=2 are supported');
    bch.m = m;
    bch.t = t;
    bch.gf = GF_init(m);
    bch.N_full = bch.gf.n;   % 2^m - 1

    % Generator polynomial: product of minimal polys of alpha^1,3,..,2t-1.
    seen = [];
    g_poly = 1;
    for i = 1:(2*t)
        coset = local_cyclotomic_coset(bch.N_full, i);
        rep = min(coset);
        if ismember(rep, seen)
            continue;
        end
        seen(end+1) = rep; %#ok<AGROW>
        m_poly = 1;
        for si = 1:numel(coset)
            s = coset(si);
            root = bch.gf.EXP(s + 1);
            m_poly = gf_polymul(bch.gf, m_poly, [root, 1]);
        end
        m_poly = bitand(round(m_poly), 1);
        g_poly = gf_polymul(bch.gf, g_poly, m_poly);
        g_poly = bitand(round(g_poly), 1);
    end
    bch.g_poly = g_poly;
    bch.n_parity = numel(g_poly) - 1;
    bch.K_full = bch.N_full - bch.n_parity;
    assert(bch.n_parity == m * t, 'bch_init: deg(g) != m*t');

    if nargin < 3 || isempty(n_short), n_short = bch.N_full; end
    if nargin < 4 || isempty(k_short), k_short = bch.K_full; end
    bch.n = n_short;
    bch.k = k_short;
    bch.shorten = bch.K_full - bch.k;
    assert(bch.shorten >= 0 && bch.N_full - bch.shorten == bch.n, ...
        'bch_init: bad shortening');
    assert(bch.n - bch.k == bch.n_parity, 'bch_init: parity count mismatch');

    bch.G = local_build_G(bch.gf, bch.N_full, bch.K_full, bch.n_parity, bch.g_poly);
    if t == 2
        bch.lut_A = build_lut_A(bch.gf);
    end
end

% ------------------------------------------------------------------------
function coset = local_cyclotomic_coset(n_full, i)
    coset = [];
    j = i;
    while ~ismember(j, coset)
        coset(end+1) = j; %#ok<AGROW>
        j = mod(j * 2, n_full);
    end
end

% ------------------------------------------------------------------------
function G = local_build_G(gf, N, K, n_minus_k, g)
    G = zeros(K, N);
    for i = 0:K-1
        dividend = zeros(1, N);
        dividend(n_minus_k + i + 1) = 1;
        [~, rem] = gf_polydivmod(gf, dividend, g);
        parity = zeros(1, n_minus_k);
        parity(1:numel(rem)) = rem;
        G(i+1, 1:n_minus_k) = parity;
        G(i+1, n_minus_k+i+1) = 1;
    end
end
