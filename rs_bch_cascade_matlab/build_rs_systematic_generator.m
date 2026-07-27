function [G, Theta_c] = build_rs_systematic_generator(gf, Theta, k_prime, n)
%BUILD_RS_SYSTEMATIC_GENERATOR  RS systematic generator via Lagrange interp.
%
%   ★ CORE INNOVATION (port of src/llosd.py:build_rs_systematic_generator).
%
%   Given the MRIP set Theta (k' positions, 1-based indices in 1..n), build a
%   k' x n RS systematic generator matrix G over GF(2^m) whose columns indexed
%   by Theta form the k' x k' identity. Each row is the Lagrange polynomial
%   T_{j_i}(x) that is 1 at locator alpha^{j_i} and 0 at the other Theta
%   locators. This is done DIRECTLY from EXP/LOG tables — NO Gaussian
%   elimination — which is exactly the low-latency claim of the paper.
%
%   NOTE: Theta uses 1-based position indices (MATLAB natural). Locator of
%   position p (1-based) is loc(p) = alpha^{p-1}, matching the Python
%   loc = EXP[0..n-1] with 0-based positions.

    assert(numel(Theta) == k_prime);
    Theta = Theta(:).';

    Theta_mask = false(1, n);
    Theta_mask(Theta) = true;
    Theta_c = find(~Theta_mask);      % 1-based parity positions, size n-k'

    loc = gf.EXP(1:n);                 % loc(p) = alpha^{p-1}

    G = zeros(k_prime, n);
    for row_idx = 1:k_prime
        G(row_idx, Theta(row_idx)) = 1;   % identity on Theta columns
    end

    if isempty(Theta_c)
        return;
    end

    a_i = loc(Theta);                 % (1, k')
    a_c = loc(Theta_c);               % (1, n-k')

    % Denominators: denom_log(i) = sum_{j'!=i} log(a_i(i) XOR a_i(j')).
    XOR_ii = bitxor(a_i.', a_i);      % (k', k') pairwise XOR (broadcast)
    kp = k_prime;
    XOR_ii(1:(kp+1):kp*kp) = 1;       % set diagonal to 1 -> log = 0
    LOG_ii = gf.LOG(XOR_ii + 1);
    denom_log = mod(sum(LOG_ii, 2), gf.n);   % (k', 1)

    % Numerators for (c, i): prod_{j'!=i} (a_c(c) XOR a_i(j')).
    XOR_ci = bitxor(a_c.', a_i);      % (n-k', k')
    LOG_ci = gf.LOG(XOR_ci + 1);      % (n-k', k')
    total_log_ci = sum(LOG_ci, 2);    % (n-k', 1)
    num_log = mod(total_log_ci - LOG_ci, gf.n);          % (n-k', k')
    exp_index = mod(num_log - denom_log.', gf.n);        % (n-k', k')
    vals = gf.EXP(exp_index + 1);                        % (n-k', k')

    % G(row i, col Theta_c(c)) = vals(c, i)  ->  G(:, Theta_c) = vals.'
    G(:, Theta_c) = vals.';
end
