function k = bch_dimension(g, n)
%BCH_DIMENSION Dimension of a binary cyclic code from a low-order polynomial.
degree = numel(g) - 1;
while degree > 0 && g(degree + 1) == 0
    degree = degree - 1;
end
k = n - degree;
end
