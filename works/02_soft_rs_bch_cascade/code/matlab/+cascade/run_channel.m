function llr = run_channel(bits, ebn0_db, rate)
%RUN_CHANNEL PAM4 modulation -> AWGN -> LLR, preserving an odd input length.
bits = double(bits(:).');
originalLength = numel(bits);
if mod(originalLength, 2) ~= 0
    bits = [bits, 0]; %#ok<AGROW>
end
sigma = cascade.sigma_from_ebn0_pam4(ebn0_db, rate);
x = cascade.bits_to_pam4(bits);
y = cascade.awgn_channel(x, sigma);
allLlr = cascade.pam4_bit_llr(y, sigma);
llr = allLlr(1:originalLength);
end
