% Elementary GF operations using EXP/LOG tables.

function y = gf_add(a, b)
    y = bitxor(a, b);
end
