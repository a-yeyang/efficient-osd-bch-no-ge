function symbols = bits_to_symbols(bits, m)
%BITS_TO_SYMBOLS Convert an LSB-first stream to GF integer symbols.
bits = double(bits(:).');
assert(mod(numel(bits), m) == 0, 'cascade:Bits:Length', 'Bit length must be a multiple of m.');
symbols = zeros(1, numel(bits) / m);
for ii = 1:numel(symbols)
    value = uint32(0);
    for bit = 0:m-1
        if bits((ii - 1) * m + bit + 1) ~= 0
            value = bitor(value, bitshift(uint32(1), bit));
        end
    end
    symbols(ii) = double(value);
end
end
