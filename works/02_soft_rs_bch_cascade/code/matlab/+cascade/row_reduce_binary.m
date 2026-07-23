function reduced = row_reduce_binary(matrix)
%ROW_REDUCE_BINARY RREF over GF(2), dropping all-zero rows.
reduced = mod(double(matrix), 2);
[nRows, nCols] = size(reduced);
row = 1;
for col = 1:nCols
    if row > nRows
        break;
    end
    pivot = find(reduced(row:end, col), 1, 'first');
    if isempty(pivot)
        continue;
    end
    pivot = pivot + row - 1;
    if pivot ~= row
        reduced([row, pivot], :) = reduced([pivot, row], :);
    end
    for rr = 1:nRows
        if rr ~= row && reduced(rr, col)
            reduced(rr, :) = bitxor(uint8(reduced(rr, :)), uint8(reduced(row, :)));
        end
    end
    row = row + 1;
end
reduced = double(reduced(any(reduced, 2), :));
end
