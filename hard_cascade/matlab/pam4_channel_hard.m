function hard = pam4_channel_hard(bits, ebn0_db, rate)
    % PAM4-modulate bits -> AWGN -> hard-decision demodulate.
    % Pads one zero bit if the input length is odd, then truncates it back off.
    % Port of hc_src/cascade_v3.py run_pam4_channel_hard.
    sigma = pam4_sigma_from_ebn0(ebn0_db, rate);
    bits = double(bits(:).');
    n_orig = numel(bits);
    if mod(n_orig, 2) ~= 0
        bits_padded = [bits, 0];
    else
        bits_padded = bits;
    end
    x = pam4_bits_to_symbols(bits_padded);
    y = awgn_channel(x, sigma);
    hard = pam4_symbols_to_bits_hard(y);
    hard = hard(1:n_orig);
end
