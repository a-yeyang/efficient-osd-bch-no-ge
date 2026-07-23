function received = awgn_channel(transmitted, sigma)
%AWGN_CHANNEL Add real iid Gaussian noise.  RNG is controlled by caller.
received = double(transmitted) + double(sigma) .* randn(size(transmitted));
end
