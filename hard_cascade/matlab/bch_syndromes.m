function [S1, S3] = bch_syndromes(gf, r)
    % Compute odd syndromes S_1, S_3 for binary BCH.
    n = gf.n;
    S1 = 0; S3 = 0;
    for j = 0:n-1
        if r(j+1) == 1
            S1 = bitxor(S1, gf.EXP(j+1));
            S3 = bitxor(S3, gf.EXP(mod(3*j, n)+1));
        end
    end
end
