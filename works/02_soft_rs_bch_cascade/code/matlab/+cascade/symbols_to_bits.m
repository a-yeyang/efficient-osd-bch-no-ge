function bits = symbols_to_bits(symbols, m)
%SYMBOLS_TO_BITS Convert GF symbols to m LSB-first binary bits each.
symbols = double(symbols(:).');
bits = zeros(1, numel(symbols) * m);
for ii = 1:numel(symbols)
    for bit = 0:m-1
        bits((ii - 1) * m + bit + 1) = bitget(uint32(symbols(ii)), bit + 1);
    end
end
end
