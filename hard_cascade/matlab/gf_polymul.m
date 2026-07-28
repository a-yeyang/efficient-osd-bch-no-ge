function out = gf_polymul(gf, a, b)
    % Polynomial multiply over GF(2^m). a, b: little-endian coeff vectors.
    a = a(:).'; b = b(:).';
    out = zeros(1, numel(a) + numel(b) - 1);
    for i = 1:numel(a)
        if a(i) == 0, continue; end
        for j = 1:numel(b)
            if b(j) == 0, continue; end
            out(i+j-1) = bitxor(out(i+j-1), gf_mul(gf, a(i), b(j)));
        end
    end
end
