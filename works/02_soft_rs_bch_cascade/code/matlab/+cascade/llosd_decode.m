function [bestCodeword, stats] = llosd_decode(code, llr, tau, use_binary_reencoding, use_early_terminate, max_bch_candidates, cache)
%LLOSD_DECODE MATLAB implementation of order-tau LLOSD/LLOSD-B.
%   It follows Work 01's reference algorithm: reliability-sort positions,
%   use a Lagrange systematic RS supercode, filter binary parity candidates,
%   and minimize correlation distance.  Outputs are binary row vectors.
if nargin < 4 || isempty(use_binary_reencoding)
    use_binary_reencoding = false;
end
if nargin < 5 || isempty(use_early_terminate)
    use_early_terminate = true;
end
if nargin < 6 || isempty(max_bch_candidates)
    max_bch_candidates = inf;
end
if nargin < 7
    cache = [];
end
llr = double(llr(:).');
assert(numel(llr) == code.n, 'cascade:LLOSD:WordLength', 'LLR vector has wrong length.');
t0 = tic;
counters = cascade.OpCounters();
n = code.n;
kPrime = n - 2 * code.t;
d = code.d_design;

rHard = double(llr < 0);
counters.fp = counters.fp + n;
[~, perm] = sort(abs(llr), 'descend');
counters.fp = counters.fp + n * floor(log2(max(n, 2)));
theta = perm(1:kPrime);
[G, thetaComplement] = cascade.build_rs_systematic_generator(code.gf, theta, kPrime, n, cache);
counters.f2m = counters.f2m + 2 * (n * n - kPrime * kPrime + kPrime);
Gpc = G(:, thetaComplement);

u0 = rHard(theta);
active = find(u0 ~= 0);
vHat = zeros(1, numel(thetaComplement));
for ii = 1:numel(active)
    vHat = cascade.gfxor(vHat, Gpc(active(ii), :));
end
counters.f2m = counters.f2m + numel(active) * numel(thetaComplement);

absLlr = abs(llr);
template = zeros(1, n);
template(theta) = u0;
bestCodeword = [];
bestDistance = inf;
nTeps = 0;
nBch = 0;
terminated = false;

for weight = 0:tau
    if weight == 0
        supports = zeros(1, 0);
    else
        supports = nchoosek(1:kPrime, weight);
    end
    if weight == 0
        nSupportRows = 1;
    else
        nSupportRows = size(supports, 1);
    end
    for supportRow = 1:nSupportRows
        if weight == 0
            support = zeros(1, 0);
        else
            support = supports(supportRow, :);
        end
        nTeps = nTeps + 1;
        parity = vHat;
        for ii = 1:numel(support)
            parity = cascade.gfxor(parity, Gpc(support(ii), :));
        end
        if use_binary_reencoding
            counters.f2 = counters.f2 + code.m * numel(thetaComplement) * max(1, numel(support));
        else
            counters.f2m = counters.f2m + numel(thetaComplement) * max(1, numel(support));
        end
        if any(parity > 1)
            continue;
        end
        nBch = nBch + 1;
        candidate = template;
        if ~isempty(support)
            candidate(theta(support)) = mod(candidate(theta(support)) + 1, 2);
        end
        candidate(thetaComplement) = parity;
        diffMask = rHard ~= candidate;
        distance = sum(absLlr(diffMask));
        counters.fp = counters.fp + n;
        if distance < bestDistance
            bestDistance = distance;
            bestCodeword = candidate;
            if use_early_terminate
                dOmega = sum(diffMask);
                matching = absLlr(~diffMask);
                K = max(0, d - dOmega - 1);
                if K == 0
                    mlOK = distance <= 0;
                else
                    matching = sort(matching, 'ascend');
                    mlOK = distance <= sum(matching(1:K));
                end
                if mlOK
                    terminated = true;
                    counters.latency_us = toc(t0) * 1e6;
                    stats = struct('counters', counters, 'n_teps', nTeps, ...
                        'n_bch_candidates', nBch, 'terminated_early', terminated);
                    return;
                end
            end
        end
        if nBch >= max_bch_candidates
            break;
        end
    end
    if nBch >= max_bch_candidates
        break;
    end
end
if isempty(bestCodeword)
    bestCodeword = rHard;
end
counters.latency_us = toc(t0) * 1e6;
stats = struct('counters', counters, 'n_teps', nTeps, ...
    'n_bch_candidates', nBch, 'terminated_early', terminated);
end
