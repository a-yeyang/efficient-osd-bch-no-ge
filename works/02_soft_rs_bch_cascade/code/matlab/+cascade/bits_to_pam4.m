function symbols = bits_to_pam4(bits)
%BITS_TO_PAM4 Gray PAM4: 00->-3, 01->-1, 11->+1, 10->+3.
bits = double(bits(:).');
assert(mod(numel(bits), 2) == 0, 'cascade:PAM4:OddBits', 'PAM4 requires an even number of bits.');
idx = 2 * bits(1:2:end) + bits(2:2:end);
levels = [-3, -1, 3, 1]; % index is binary b1b0: 00,01,10,11
symbols = levels(idx + 1);
end
