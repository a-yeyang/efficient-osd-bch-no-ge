function [c_short, ok] = shortened_rs_decode(srs, r_short)
    % Re-insert known zeros -> full-length BM decode -> drop known zeros.
    r_full = [r_short(1:srs.n_parity), zeros(1, srs.shorten), r_short(srs.n_parity+1:end)];
    [c_dec_full, ok] = rs_bm_decode(srs.full, r_full);
    c_short = [c_dec_full(1:srs.n_parity), c_dec_full(srs.n_parity + srs.shorten + 1:end)];
end
