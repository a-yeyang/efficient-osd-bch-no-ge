function inverse = gf2_inverse(matrix)
%GF2_INVERSE Invert a nonsingular square binary matrix over GF(2).
matrix = mod(double(matrix), 2);
[nRows, nCols] = size(matrix);
assert(nRows == nCols, 'cascade:GF2:NotSquare', 'Matrix must be square.');
augmented = [uint8(matrix), uint8(eye(nRows))];
for col = 1:nRows
    pivot = find(augmented(col:end, col), 1, 'first');
    if isempty(pivot)
        error('cascade:GF2:Singular', 'Matrix is singular over GF(2).');
    end
    pivot = pivot + col - 1;
    if pivot ~= col
        augmented([col, pivot], :) = augmented([pivot, col], :);
    end
    for rr = 1:nRows
        if rr ~= col && augmented(rr, col)
            augmented(rr, :) = bitxor(augmented(rr, :), augmented(col, :));
        end
    end
end
inverse = double(augmented(:, nRows + 1:end));
end
