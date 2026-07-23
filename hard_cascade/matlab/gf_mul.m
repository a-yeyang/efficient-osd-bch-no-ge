function y = gf_mul(gf, a, b)
    % GF(2^m) multiplication via EXP/LOG tables
    if a == 0 || b == 0
        y = 0;
        return;
    end
    la = gf.LOG(a + 1);
    lb = gf.LOG(b + 1);
    y = gf.EXP(la + lb + 1);
end
