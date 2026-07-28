function msg = bch_extract_message(bch, cw_short)
    msg = cw_short(bch.n_parity + 1:end);
end
