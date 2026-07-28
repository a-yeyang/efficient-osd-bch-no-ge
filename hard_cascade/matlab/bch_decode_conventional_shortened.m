function [c_short, ok] = bch_decode_conventional_shortened(bch, r_short)
    % Shortened-code wrapper: reinsert known zeros, decode full length, drop them.
    % t=1 reduces to the closed form p = log(S1) (Conventional == Direct);
    % t=2 reuses the tested full-length bch_decode_conventional.m.
    r_full = [r_short(1:bch.n_parity), zeros(1, bch.shorten), r_short(bch.n_parity+1:end)];
    if bch.t == 1
        [S1, ~] = bch_syndromes(bch.gf, r_full);
        c_full = r_full;
        if S1 ~= 0
            p = mod(bch.gf.LOG(S1 + 1), bch.N_full);
            c_full(p+1) = 1 - c_full(p+1);
        end
        ok = true;
    else
        [c_full, ok] = bch_decode_conventional(bch.gf, r_full);
    end
    c_short = [c_full(1:bch.n_parity), c_full(bch.n_parity + bch.shorten + 1:end)];
end
