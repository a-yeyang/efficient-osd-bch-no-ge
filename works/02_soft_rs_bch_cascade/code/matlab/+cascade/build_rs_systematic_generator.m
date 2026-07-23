function [G, thetaComplement] = build_rs_systematic_generator(gf, theta, k_prime, n, cache)
%BUILD_RS_SYSTEMATIC_GENERATOR Lagrange-form systematic RS generator matrix.
%   THETA and THETACOMPLEMENT use MATLAB one-based codeword positions.  This
%   is a deliberately direct form of the Python reference so it remains
%   readable and toolbox-free; the code sizes used here have n-k' <= 6.
theta = double(theta(:).');
if nargin < 5
    cache = [];
end
assert(numel(theta) == k_prime, 'cascade:LLOSD:MRIPLength', 'Theta must contain k'' positions.');
assert(numel(unique(theta)) == k_prime && all(theta >= 1) && all(theta <= n), ...
    'cascade:LLOSD:MRIPPositions', 'Theta must be unique positions in 1..n.');
thetaMask = false(1, n);
thetaMask(theta) = true;
thetaComplement = find(~thetaMask);
G = zeros(k_prime, n);
for row = 1:k_prime
    G(row, theta(row)) = 1;
end
if isempty(thetaComplement)
    return;
end

if ~isempty(cache)
    % Scheme B's actual shared path: the pairwise locator table is built once
    % and all Lagrange numerator/denominator requests pass through CACHE.
    % This is intentionally not a bookkeeping-only call: every generator
    % matrix element below is evaluated using cached field differences and
    % cached denominator products.
    cache.build_pairwise_diff();
    for row = 1:k_prime
        for cc = 1:numel(thetaComplement)
            column = thetaComplement(cc);
            G(row, column) = cache.lagrange_basis(theta(row), theta, column);
        end
    end
    return;
end

locators = gf.EXP(1:n); % alpha^(position-1)
for row = 1:k_prime
    denominator = 1;
    aRow = locators(theta(row));
    for other = 1:k_prime
        if other ~= row
            denominator = gf.mul(denominator, cascade.gfxor(aRow, locators(theta(other))));
        end
    end
    for cc = 1:numel(thetaComplement)
        column = thetaComplement(cc);
        numerator = 1;
        for other = 1:k_prime
            if other ~= row
                numerator = gf.mul(numerator, cascade.gfxor(locators(column), locators(theta(other))));
            end
        end
        G(row, column) = gf.div(numerator, denominator);
    end
end
end
