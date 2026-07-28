function s = pam4_sigma_from_ebn0(ebn0_db, rate)
    % sigma^2 = E_s_avg * rate / (4 * Eb/N0), E_s_avg = 5.0 for PAM4.
    E_S_AVG = 5.0;
    ebn0_lin = 10.0 ^ (ebn0_db / 10.0);
    sigma2 = E_S_AVG * rate / (4.0 * ebn0_lin);
    s = sqrt(sigma2);
end
