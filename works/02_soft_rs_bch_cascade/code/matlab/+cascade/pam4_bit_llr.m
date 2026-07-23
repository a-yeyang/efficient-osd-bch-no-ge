function llr = pam4_bit_llr(received, sigma, use_max_log)
%PAM4_BIT_LLR Per-bit log P(b=0|y)/P(b=1|y), stable full log-sum-exp.
if nargin < 3
    use_max_log = false;
end
y = double(received(:).');
inv2s2 = 1 / (2 * double(sigma)^2);
dm3 = -inv2s2 .* (y + 3).^2;
dm1 = -inv2s2 .* (y + 1).^2;
dp1 = -inv2s2 .* (y - 1).^2;
dp3 = -inv2s2 .* (y - 3).^2;
if use_max_log
    b1 = max(dm3, dm1) - max(dp1, dp3);
    b0 = max(dm3, dp3) - max(dm1, dp1);
else
    b1 = cascade.logsumexp2(dm3, dm1) - cascade.logsumexp2(dp1, dp3);
    b0 = cascade.logsumexp2(dm3, dp3) - cascade.logsumexp2(dm1, dp1);
end
llr = zeros(1, 2 * numel(y));
llr(1:2:end) = b1;
llr(2:2:end) = b0;
end
