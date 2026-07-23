function [pivots, inverse_submatrix] = gf2_independent_columns(generator)
%GF2_INDEPENDENT_COLUMNS Pick k independent columns of a k-by-n GF(2) G.
original = mod(double(generator), 2);
[k, n] = size(original);
work = uint8(original);
pivots = zeros(1, k);
row = 1;
for col = 1:n
    if row > k
        break;
    end
    pivot = find(work(row:end, col), 1, 'first');
    if isempty(pivot)
        continue;
    end
    pivot = pivot + row - 1;
    if pivot ~= row
        work([row, pivot], :) = work([pivot, row], :);
    end
    for rr = 1:k
        if rr ~= row && work(rr, col)
            work(rr, :) = bitxor(work(rr, :), work(row, :));
        end
    end
    pivots(row) = col;
    row = row + 1;
end
if row <= k
    error('cascade:GF2:RankDeficient', 'Generator matrix rank is smaller than k.');
end
inverse_submatrix = cascade.gf2_inverse(original(:, pivots));
end
