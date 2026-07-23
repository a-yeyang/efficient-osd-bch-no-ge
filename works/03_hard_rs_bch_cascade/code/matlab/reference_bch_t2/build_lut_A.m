function lut_A = build_lut_A(gf)
    % Build LUT for BCH Direct root finding (Lagendijk Table I).
    % A(X) = X^2 + X + k has roots {X, X+1} when k = X^2 + X (=X^2 XOR X in char 2)
    % lut_A.roots(k+1, :) = [X, X+1] when valid; lut_A.valid(k+1) = true
    n_field = bitshift(1, gf.m);
    lut_A.roots = zeros(n_field, 2);
    lut_A.valid = false(n_field, 1);
    for X = 0:n_field-1
        X2 = gf_mul(gf, X, X);
        k = bitxor(X2, X);
        if ~lut_A.valid(k+1)
            lut_A.roots(k+1, 1) = X;
            lut_A.roots(k+1, 2) = bitxor(X, 1);
            lut_A.valid(k+1) = true;
        end
    end
end
