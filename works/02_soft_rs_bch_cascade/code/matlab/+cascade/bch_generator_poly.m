function g = bch_generator_poly(gf, t)
%BCH_GENERATOR_POLY Narrow-sense primitive binary BCH generator polynomial.
%   Coefficients are lowest power first and agree with the Python reference.
seen = false(1, gf.n);
g = 1;
for ii = 1:(2 * t)
    coset = [];
    value = mod(ii, gf.n);
    while ~ismember(value, coset)
        coset(end + 1) = value; %#ok<AGROW>
        value = mod(2 * value, gf.n);
    end
    representative = min(coset);
    if seen(representative + 1)
        continue;
    end
    seen(representative + 1) = true;
    minimal = 1;
    for jj = 1:numel(coset)
        root = gf.alpha(coset(jj));
        minimal = gf.poly_mul(minimal, [root, 1]);
    end
    if any(~ismember(minimal, [0, 1]))
        error('cascade:BCH:NonBinaryMinimalPolynomial', ...
            'The selected primitive polynomial did not generate a binary BCH minimal polynomial.');
    end
    g = gf.poly_mul(g, minimal);
    g = mod(g, 2);
end
end
