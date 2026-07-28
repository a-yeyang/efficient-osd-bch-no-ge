function out = pam4_bits_to_symbols(bits)
    % Map a binary vector (length 2N) to N PAM4 symbols.
    % Gray mapping (b1=MSB, b0=LSB): 00->-3, 01->-1, 11->+1, 10->+3.
    bits = double(bits(:).');
    assert(mod(numel(bits), 2) == 0, 'pam4_bits_to_symbols: odd bit length');
    n_syms = numel(bits) / 2;
    out = zeros(1, n_syms);
    for i = 1:n_syms
        b1 = bits(2*i - 1);
        b0 = bits(2*i);
        if b1 == 0 && b0 == 0
            out(i) = -3.0;
        elseif b1 == 0 && b0 == 1
            out(i) = -1.0;
        elseif b1 == 1 && b0 == 1
            out(i) = +1.0;
        else
            out(i) = +3.0;
        end
    end
end
