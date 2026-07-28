function y = awgn_channel(x, sigma)
    % Add real AWGN noise with std sigma. Uses MATLAB's global rng state,
    % so callers control reproducibility via rng(seed) before invoking.
    y = x + sigma * randn(size(x));
end
