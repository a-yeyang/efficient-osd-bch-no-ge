function [best_c, info] = llosd_decode(code, L, tau, counters)
%LLOSD_DECODE  Order-tau Low-Latency OSD for BCH (port of src/llosd.py).
%
%   ★ THE PAPER'S METHOD. Because BCH is a binary subcode of RS, we build the
%   RS (n,k') systematic generator G_RS via Lagrange interpolation (NO Gaussian
%   elimination), enumerate test-error-patterns (TEPs) on the k' most-reliable
%   positions Theta, re-encode each as an RS codeword, and KEEP ONLY those whose
%   parity symbols are binary (Theorem 2) — those lie in the BCH subcode. Among
%   the surviving binary candidates, pick minimum correlation distance, with the
%   eq.(14) ML early-termination.
%
%   code : BCHCode object
%   L    : LLR vector (length n)
%   tau  : decoding order
%   counters : OpCounters handle (accumulated in place)
%
%   Returns best_c (length-n 0/1 codeword) and info struct.

    if nargin < 4 || isempty(counters), counters = OpCounters(); end
    n = code.n; m = code.m; t = code.t;
    d = code.d_design; gf = code.gf;
    k_prime = n - 2*t;
    L = double(L(:).');

    t0 = tic;

    r_hard = double(L < 0);
    counters.fp = counters.fp + n;

    % Sort positions by |L| descending -> Theta = k' most reliable (1-based).
    absL = abs(L);
    [~, perm] = sort(absL, 'descend');
    counters.fp = counters.fp + n * floor(log2(max(n, 2)));
    Theta = perm(1:k_prime);

    % ★ Build RS systematic generator via Lagrange (no GE).
    [G_RS, Theta_c] = build_rs_systematic_generator(gf, Theta, k_prime, n);
    counters.f2m = counters.f2m + 2 * (n*n - k_prime*k_prime + k_prime);

    G_pc = G_RS(:, Theta_c);                 % (k', n-k') parity slice
    ncol = numel(Theta_c);

    u0 = r_hard(Theta);                      % binary MRIP values
    active = find(u0 ~= 0);
    if ~isempty(active)
        v_hat0_pc = G_pc(active(1), :);
        for q = 2:numel(active)
            v_hat0_pc = bitxor(v_hat0_pc, G_pc(active(q), :));
        end
    else
        v_hat0_pc = zeros(1, ncol);
    end
    counters.f2m = counters.f2m + numel(active) * ncol;

    best_c = [];
    best_D = inf;
    n_teps = 0;
    n_bch = 0;
    terminated_early = false;

    for w = 0:tau
        if w == 0
            supports = zeros(1, 0);           % single empty support
        else
            supports = nchoosek(1:k_prime, w); % each row a support (1-based)
        end
        for si = 1:size(supports, 1)
            n_teps = n_teps + 1;
            support = supports(si, :);
            if isempty(support)
                v_parity = v_hat0_pc;
            else
                v_parity = v_hat0_pc;
                for q = 1:numel(support)
                    v_parity = bitxor(v_parity, G_pc(support(q), :));
                end
            end
            % Per-TEP re-encoding cost is charged to f2 at return (binary
            % re-encoding, the paper's default), based on n_teps processed.

            % Binary filter (Theorem 2): any parity symbol > 1 => non-binary.
            if any(v_parity > 1)
                continue;
            end
            n_bch = n_bch + 1;

            c_hat = zeros(1, n);
            c_hat(Theta) = u0;
            if ~isempty(support)
                c_hat(Theta(support)) = bitxor(c_hat(Theta(support)), 1);
            end
            c_hat(Theta_c) = v_parity;

            diff_mask = (r_hard ~= c_hat);
            D = sum(absL(diff_mask));
            counters.fp = counters.fp + n;
            if D < best_D
                best_D = D;
                best_c = c_hat;
                % eq.(14) ML early termination.
                d_omega = sum(diff_mask);
                match_mask = ~diff_mask;
                if any(match_mask)
                    absL_match_sorted = sort(absL(match_mask), 'ascend');
                    K = max(0, d - d_omega - 1);
                    if K == 0
                        ml_ok = (D <= 0);
                    else
                        K = min(K, numel(absL_match_sorted));
                        ml_ok = (D <= sum(absL_match_sorted(1:K)));
                    end
                    if ml_ok
                        terminated_early = true;
                        % Per-TEP binary re-encoding cost -> f2 (Python default
                        % use_binary_reencoding=True: f2 += n_teps*(n-k')*m).
                        counters.f2 = counters.f2 + n_teps * (n - k_prime) * m;
                        counters.latency_us = toc(t0) * 1e6;
                        best_c = double(best_c);
                        info = struct('counters', counters, 'n_teps', n_teps, ...
                            'n_bch_candidates', n_bch, 'terminated_early', true);
                        return;
                    end
                end
            end
        end
    end

    if isempty(best_c)
        best_c = r_hard;
    end
    % Per-TEP binary re-encoding cost -> f2 (Python default
    % use_binary_reencoding=True: f2 += n_teps*(n-k')*m).
    counters.f2 = counters.f2 + n_teps * (n - k_prime) * m;
    counters.latency_us = toc(t0) * 1e6;
    best_c = double(best_c);
    info = struct('counters', counters, 'n_teps', n_teps, ...
        'n_bch_candidates', n_bch, 'terminated_early', false);
end
