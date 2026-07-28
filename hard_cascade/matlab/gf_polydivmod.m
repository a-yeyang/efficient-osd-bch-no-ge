function [quot, rem] = gf_polydivmod(gf, num, den)
    % Polynomial division over GF(2^m). num, den: little-endian coeff vectors.
    num = num(:).'; den = den(:).';
    while ~isempty(den) && den(end) == 0
        den(end) = [];
    end
    if isempty(den)
        error('gf_polydivmod: zero divisor');
    end
    lead_inv = gf_div(gf, 1, den(end));
    ld = numel(den);
    qlen = max(0, numel(num) - ld + 1);
    quot = zeros(1, qlen);
    rem = num;
    for i = qlen-1:-1:0
        if numel(rem) <= i + ld - 1
            continue;
        end
        coeff = gf_mul(gf, rem(i + ld - 1 + 1), lead_inv);
        quot(i+1) = coeff;
        for j = 0:ld-1
            rem(i + j + 1) = bitxor(rem(i + j + 1), gf_mul(gf, coeff, den(j+1)));
        end
    end
    while ~isempty(rem) && rem(end) == 0
        rem(end) = [];
    end
end
