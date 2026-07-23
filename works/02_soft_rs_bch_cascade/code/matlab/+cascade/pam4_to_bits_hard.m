function bits = pam4_to_bits_hard(received)
%PAM4_TO_BITS_HARD Nearest-neighbour hard decision for the Gray constellation.
received = double(received(:).');
levels = [-3, -1, 1, 3];
pairs = [0, 0; 0, 1; 1, 1; 1, 0];
bits = zeros(1, 2 * numel(received));
for ii = 1:numel(received)
    [~, pos] = min(abs(levels - received(ii)));
    bits(2 * ii - 1:2 * ii) = pairs(pos, :);
end
end
