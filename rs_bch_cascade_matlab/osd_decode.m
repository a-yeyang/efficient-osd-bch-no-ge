function [c_hat, info] = osd_decode(code, L, tau, counters)
%OSD_DECODE  Traditional order-tau OSD for BCH (port of src/osd.py).
%
%   ★ CONTRAST GROUP for the innovation. Same performance as LLOSD, but builds
%   the systematic generator by GAUSSIAN ELIMINATION on the |L|-sorted G, rather
%   than by Lagrange interpolation. Used to demonstrate that Lagrange lowers the
%   operation count / latency at equal BER.
%
%   Pipeline: sort columns of G by |L| desc -> GE to expose k x k identity in
%   the first k columns -> initial message = first k sorted hard bits ->
%   enumerate TEPs weight<=tau -> re-encode, correlation distance, keep min,
%   eq.(14) ML early terminate -> un-permute.

    if nargin < 4 || isempty(counters), counters = OpCounters(); end
    n = code.n; k = code.k; d = code.d_design;
    L = double(L(:).');

    t0 = tic;

    r_hard = double(L < 0);
    counters.fp = counters.fp + n;

    absL = abs(L);
    [~, perm] = sort(absL, 'descend');       % 1-based
    counters.fp = counters.fp + n * floor(log2(max(n, 2)));

    Gp = code.G(:, perm);                     % permuted generator
    [G_sys, col_perm, ge_f2] = gaussianElimBinary(Gp);
    counters.f2 = counters.f2 + ge_f2;
    % G_sys is already in effective column order; perm_eff = perm(col_perm).
    perm_eff = perm(col_perm);

    r_sorted = r_hard(perm_eff);
    L_sorted = L(perm_eff);
    absL_sorted = abs(L_sorted);

    f_init = r_sorted(1:k);

    best_c_perm = [];
    best_D = inf;
    n_teps = 0;
    terminated_early = false;

    for w = 0:tau
        if w == 0
            supports = zeros(1, 0);
        else
            supports = nchoosek(1:k, w);
        end
        for si = 1:size(supports, 1)
            n_teps = n_teps + 1;
            support = supports(si, :);
            f_omega = f_init;
            if ~isempty(support)
                f_omega(support) = bitxor(f_omega(support), 1);
            end
            c_sorted = mod(f_omega * G_sys, 2);
            counters.f2 = counters.f2 + k * n;

            diff = (r_sorted ~= c_sorted);
            D = sum(absL_sorted(diff));
            counters.fp = counters.fp + n;
            if D < best_D
                best_D = D;
                best_c_perm = c_sorted;
                d_omega = sum(diff);
                match = ~diff;
                if any(match)
                    absL_match_sorted = sort(absL_sorted(match), 'ascend');
                    K = max(0, d - d_omega - 1);
                    if K == 0
                        ml_ok = (D <= 0);
                    else
                        K = min(K, numel(absL_match_sorted));
                        ml_ok = (D <= sum(absL_match_sorted(1:K)));
                    end
                    if ml_ok
                        terminated_early = true;
                        break;
                    end
                end
            end
        end
        if terminated_early, break; end
    end

    c_hat = zeros(1, n);
    c_hat(perm_eff) = best_c_perm;
    counters.latency_us = toc(t0) * 1e6;
    info = struct('counters', counters, 'n_teps', n_teps, ...
        'n_bch_candidates', n_teps);
end


function [G_sys, col_perm, f2_ops] = gaussianElimBinary(G_perm)
%GAUSSIANELIMBINARY  Row-reduce permuted k x n binary G so first k cols = I.
    G = mod(G_perm, 2);
    [k, n] = size(G);
    col_perm = 1:n;
    f2_ops = 0;
    for i = 0:k-1
        piv_col = [];
        for c = i:n-1
            if G(i+1, c+1) ~= 0, piv_col = c; break; end
        end
        if isempty(piv_col)
            piv_row = []; pc = [];
            for r = i+1:k-1
                for c = i:n-1
                    if G(r+1, c+1) ~= 0, piv_row = r; pc = c; break; end
                end
                if ~isempty(piv_row), break; end
            end
            if isempty(piv_row)
                error('osd:degenerate', 'degenerate G during GE');
            end
            G([i+1, piv_row+1], :) = G([piv_row+1, i+1], :);
            piv_col = pc;
        end
        if piv_col ~= i
            G(:, [i+1, piv_col+1]) = G(:, [piv_col+1, i+1]);
            col_perm([i+1, piv_col+1]) = col_perm([piv_col+1, i+1]);
        end
        for r = 0:k-1
            if r ~= i && G(r+1, i+1) ~= 0
                G(r+1, :) = bitxor(G(r+1, :), G(i+1, :));
                f2_ops = f2_ops + n;
            end
        end
    end
    G_sys = G;
end
