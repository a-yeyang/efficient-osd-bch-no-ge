function sigma = sigma_from_ebn0_pam4(ebn0_db, rate)
%SIGMA_FROM_EBN0_PAM4 Python-reference convention: Es=5, N0=2*sigma^2.
ebn0 = 10.^(double(ebn0_db) / 10);
sigma = sqrt(5 .* double(rate) ./ (4 .* ebn0));
end
