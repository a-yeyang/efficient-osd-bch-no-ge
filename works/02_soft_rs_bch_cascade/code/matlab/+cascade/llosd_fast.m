function [codeword, stats] = llosd_fast(code, llr, tau, use_binary_reencoding, use_early_terminate, cache)
%LLOSD_FAST API-compatible MATLAB fallback for the NumPy/Numba fast wrapper.
%   R2022b Base has no JIT dependency here; the vector/small-parity MATLAB
%   implementation is deterministic and implements the same candidate logic.
if nargin < 4
    use_binary_reencoding = true;
end
if nargin < 5
    use_early_terminate = true;
end
if nargin < 6
    cache = [];
end
[codeword, stats] = cascade.llosd_decode(code, llr, tau, ...
    use_binary_reencoding, use_early_terminate, inf, cache);
end
