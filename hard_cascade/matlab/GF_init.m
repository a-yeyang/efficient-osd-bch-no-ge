% GF(2^m) finite field arithmetic — EXP/LOG tables.
%
% Usage:
%   gf = GF_init(m, prim_poly)
%   gf = GF_init(8)           % uses default primitive polynomial
%
% All ops are integer arithmetic on elements 0..2^m-1.

function gf = GF_init(m, prim_poly)
    if nargin < 2
        prims = containers.Map({1,2,3,4,5,6,7,8}, ...
            {bin2dec('11'), bin2dec('111'), bin2dec('1011'), ...
             bin2dec('10011'), bin2dec('100101'), bin2dec('1000011'), ...
             bin2dec('10001001'), bin2dec('100011101')});
        prim_poly = prims(m);
    end
    gf.m = m;
    gf.n = bitshift(1, m) - 1;   % 2^m - 1
    gf.prim = prim_poly;
    gf.EXP = zeros(1, 2 * gf.n + 2);
    gf.LOG = -ones(1, bitshift(1, m));
    x = 1;
    for i = 0:gf.n-1
        gf.EXP(i + 1) = x;
        gf.LOG(x + 1) = i;
        x = bitshift(x, 1);
        if bitand(x, bitshift(1, m)) ~= 0
            x = bitxor(x, prim_poly);
        end
    end
    for i = gf.n:2*gf.n+1
        gf.EXP(i + 1) = gf.EXP(i - gf.n + 1);
    end
end
