function [lambda, L] = berlekamp_massey(gf, syndromes)
%BERLEKAMP_MASSEY Error-locator polynomial for S_1..S_N over GF(2^m).
syndromes = double(syndromes(:).');
L = 0;
lambda = 1;
B = 1;
b = 1;
shift = 1;
for kk = 1:numel(syndromes)
    delta = syndromes(kk);
    for ii = 1:L
        if ii + 1 <= numel(lambda) && lambda(ii + 1) ~= 0
            delta = cascade.gfxor(delta, gf.mul(lambda(ii + 1), syndromes(kk - ii)));
        end
    end
    if delta == 0
        shift = shift + 1;
        continue;
    end
    coefficient = gf.div(delta, b);
    shiftedB = [zeros(1, shift), B];
    newLength = max(numel(lambda), numel(shiftedB));
    candidate = [lambda, zeros(1, newLength - numel(lambda))];
    shiftedB = [shiftedB, zeros(1, newLength - numel(shiftedB))];
    for ii = 1:newLength
        candidate(ii) = cascade.gfxor(candidate(ii), gf.mul(coefficient, shiftedB(ii)));
    end
    if 2 * L <= kk - 1
        Lnew = kk - L;
        B = lambda;
        b = delta;
        lambda = candidate;
        L = Lnew;
        shift = 1;
    else
        lambda = candidate;
        shift = shift + 1;
    end
end
while numel(lambda) > 1 && lambda(end) == 0
    lambda(end) = [];
end
end
