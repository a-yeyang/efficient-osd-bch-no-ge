function [message_hat, result] = cascade_scheme_b_decode(codec, llr, cache, counters)
%CASCADE_SCHEME_B_DECODE Scheme A decoding plus reusable Lagrange accounting.
if nargin < 3 || isempty(cache)
    cache = cascade.LagrangeCache(codec.rs.gf, codec.rs.n);
end
if nargin < 4 || isempty(counters)
    counters = cascade.OpCounters();
end
[message_hat, result] = codec.decode_scheme_a(llr, counters, cache);
% Use observable cache reuse rather than a fixed nominal reduction.  The
% inner Lagrange generator calls `lagrange_basis` for every parity entry;
% denominator hits are actual cache reuses, while RS BM/Chien obtains its
% alpha locators through the same cache.  `ops_saved` is in the same canonical
% accounting unit used by LLOSD's generator-cost formula.
cacheStats = cache.stats();
savings = min(counters.f2m, cacheStats.ops_saved);
counters.f2m = counters.f2m - savings;
result.counters = counters;
result.cache_stats = cache.stats();
result.cache_accounted_savings = savings;
end
