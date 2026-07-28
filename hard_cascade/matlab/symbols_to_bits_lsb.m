function bits = symbols_to_bits_lsb(symbols, m)
    % Serialize GF(2^m) symbols to bits, LSB first (matches cascade_encode/decode).
    n = numel(symbols);
    bits = zeros(1, n*m);
    for i = 0:n-1
        si = symbols(i+1);
        for b = 0:m-1
            bits(i*m + b + 1) = bitand(bitshift(si, -b), 1);
        end
    end
end
