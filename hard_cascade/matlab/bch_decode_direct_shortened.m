function [c_short, ok] = bch_decode_direct_shortened(bch, r_short)
    % Shortened-code wrapper for the Direct decoder.
    % t=1: no quadratic/no LUT — closed form p = log(S1) (identical to Conventional).
    % t=2: reuses the tested full-length bch_decode_direct.m + bch.lut_A.
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
        [c_full, ok] = bch_decode_direct(bch.gf, r_full, bch.lut_A);
    end
    c_short = [c_full(1:bch.n_parity), c_full(bch.n_parity + bch.shorten + 1:end)];
end
