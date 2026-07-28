function msg = rs_extract_message(rs, c)
    msg = c(rs.n - rs.k + 1 : end);
end
