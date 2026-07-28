function out = pam4_symbols_to_bits_hard(y)
    % Hard-decision demodulator: nearest of {-3,-1,+1,+3} -> Gray bit pair.
    y = y(:).';
    out = zeros(1, 2*numel(y));
    levels = [-3.0, -1.0, +1.0, +3.0];
    for i = 1:numel(y)
        [~, mi] = min(abs(levels - y(i)));
        level = levels(mi);
        switch level
            case -3.0, b1 = 0; b0 = 0;
            case -1.0, b1 = 0; b0 = 1;
            case +1.0, b1 = 1; b0 = 1;
            otherwise, b1 = 1; b0 = 0;   % +3
        end
        out(2*i - 1) = b1;
        out(2*i)     = b0;
    end
end
