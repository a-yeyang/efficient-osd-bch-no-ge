function [solution, ok] = gf_solve(gf, matrix, rhs)
%GF_SOLVE Solve A*x=b over GF(2^m) by Gauss-Jordan elimination.
matrix = double(matrix);
rhs = double(rhs(:));
[nRows, nCols] = size(matrix);
if nRows ~= nCols || numel(rhs) ~= nRows
    error('cascade:GF:BadLinearSystem', 'A must be square and conformable with b.');
end
augmented = [matrix, rhs];
ok = true;
for col = 1:nCols
    pivot = find(augmented(col:end, col) ~= 0, 1, 'first');
    if isempty(pivot)
        solution = zeros(nCols, 1);
        ok = false;
        return;
    end
    pivot = pivot + col - 1;
    if pivot ~= col
        augmented([col, pivot], :) = augmented([pivot, col], :);
    end
    pivotInverse = gf.inv(augmented(col, col));
    for jj = col:(nCols + 1)
        augmented(col, jj) = gf.mul(augmented(col, jj), pivotInverse);
    end
    for rr = 1:nRows
        if rr == col || augmented(rr, col) == 0
            continue;
        end
        factor = augmented(rr, col);
        for jj = col:(nCols + 1)
            augmented(rr, jj) = cascade.gfxor(augmented(rr, jj), gf.mul(factor, augmented(col, jj)));
        end
    end
end
solution = augmented(:, end);
end
