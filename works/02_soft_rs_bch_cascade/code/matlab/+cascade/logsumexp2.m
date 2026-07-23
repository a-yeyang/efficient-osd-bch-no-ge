function out = logsumexp2(a, b)
%LOGSUMEXP2 Stable log(exp(a)+exp(b)) for matching-size arrays.
m = max(a, b);
out = m + log1p(exp(-abs(a - b)));
end
