function y = gf_div(gf, a, b)
    % GF(2^m) division
    if a == 0
        y = 0; return;
    end
    if b == 0
        error('divide by zero in GF');
    end
    la = gf.LOG(a + 1);
    lb = gf.LOG(b + 1);
    y = gf.EXP(mod(la - lb + gf.n, gf.n) + 1);
end
