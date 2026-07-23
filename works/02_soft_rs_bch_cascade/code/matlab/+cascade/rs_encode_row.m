function codeword = rs_encode_row(gf, message, generator)
%RS_ENCODE_ROW Encode binary message coefficients through an RS generator.
message = double(message(:).');
[k, n] = size(generator);
assert(numel(message) == k, 'cascade:RS:MessageLength', 'Message length does not match generator.');
codeword = zeros(1, n);
for ii = 1:k
    if message(ii) ~= 0
        codeword = cascade.gfxor(codeword, generator(ii, :));
    end
end
end
