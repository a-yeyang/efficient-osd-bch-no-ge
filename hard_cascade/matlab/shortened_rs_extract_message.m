function msg = shortened_rs_extract_message(srs, c_short)
    msg = c_short(srs.n_parity + 1:end);
end
