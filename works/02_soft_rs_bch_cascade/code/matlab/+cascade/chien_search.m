function positions = chien_search(gf, lambda, n, cache)
%CHIEN_SEARCH Return one-based word positions p whose alpha^(-p0) is a root.
lambda = double(lambda(:).');
if nargin < 4
    cache = [];
end
positions = [];
for position = 1:n
    p0 = position - 1;
    value = 0;
    for degree = 0:numel(lambda)-1
        if lambda(degree + 1) ~= 0
            exponentPosition = mod((n - p0) * degree, gf.n) + 1;
            if isempty(cache)
                locator = gf.alpha(exponentPosition - 1);
            else
                locator = cache.alpha_at(exponentPosition);
            end
            value = cascade.gfxor(value, gf.mul(lambda(degree + 1), locator));
        end
    end
    if value == 0
        positions(end + 1) = position; %#ok<AGROW>
    end
end
end
