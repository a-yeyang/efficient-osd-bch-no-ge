function value = theoretical_puncture(code, nRandomPuncturings, seed)
%THEORETICAL_PUNCTURE Direct small-code analogue of fig02's Python helper.
% This intentionally returns [] for k>22, matching the Python implementation.

if nargin < 2, nRandomPuncturings = 10000; end
if nargin < 3, seed = 0; end
if code.k > 22, value = []; return; end
rng(seed,'twister');
nWords = 2^code.k;
msgs = zeros(nWords,code.k,'uint8');
ids = uint32(0:nWords-1).';
for j = 1:code.k, msgs(:,j) = uint8(bitget(ids,j)); end
cws = mod(double(msgs)*code.G,2);
kPrime = code.n-2*code.t;
counts = zeros(1,nRandomPuncturings);
for q = 1:nRandomPuncturings
    punctured = randperm(code.n,code.n-kPrime);
    keep = true(1,code.n); keep(punctured) = false;
    counts(q) = nnz(sum(cws(:,keep),2) <= code.t);
end
value = mean(counts);
end
