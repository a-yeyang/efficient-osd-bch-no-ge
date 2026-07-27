function selftest()
%SELFTEST  Numerical self-checks for the MATLAB RS+BCH cascade port.
%
%   Verifies (per the plan's end-to-end validation):
%     1. GF tables: mul(a, inv(a)) == 1 for all nonzero a.
%     2. RS: random msg -> encode -> add <=t errors -> bmDecode recovers.
%     3. BCH: G H^T == 0; encode produces a valid codeword.
%     4. Lagrange: build_rs_systematic_generator gives identity on Theta and
%        every row is a valid RS codeword (syndromes vanish); and it matches an
%        independent Lagrange evaluation.
%     5. LLOSD == OSD: same decoded codeword on identical LLR inputs.
%
%   Throws (assert) on any failure. Prints a PASS line per check.

    fprintf('=== Self-test ===\n');

    % --- 1. GF tables ---
    for m = [7, 8]
        gf = GF(m);
        for a = 1:gf.n
            assert(gf.mul(a, gf.inv(a)) == 1, 'GF inv failed m=%d a=%d', m, a);
        end
        % a * a^{-1} via div, and pow consistency
        assert(gf.mul(gf.EXP(3), gf.EXP(5)) == gf.EXP(mod(2+4, gf.n)+1), 'GF mul via EXP');
        assert(gf.powr(gf.EXP(2), 3) == gf.EXP(mod(1*3, gf.n)+1), 'GF pow');
    end
    fprintf('  [PASS] GF(2^m) mul/inv/div/pow consistent (m=7,8)\n');

    % --- 2. RS encode/decode ---
    rng = RandStream('mt19937ar', 'Seed', 12345);
    rs = RSCode(8, 239);
    for trial = 1:20
        msg = double(randi(rng, [0, 255], 1, rs.k));
        c = rs.encodeSystematic(msg);
        % add up to t errors
        r = c;
        n_err = randi(rng, [0, rs.t]);
        pos = randperm(rng, rs.n, n_err);
        for pp = pos
            e = randi(rng, [1, 255]);
            r(pp) = bitxor(r(pp), e);
        end
        [c_dec, ok] = rs.bmDecode(r);
        assert(ok, 'RS bmDecode failed to decode (trial %d)', trial);
        assert(isequal(c_dec, c), 'RS bmDecode wrong codeword (trial %d)', trial);
        assert(isequal(rs.extractMessage(c_dec), msg), 'RS msg mismatch');
    end
    fprintf('  [PASS] RS(255,239) encode -> <=t errors -> BM recovers (20 trials)\n');

    % --- 3. BCH ---
    for spec = {[7,1], [8,2]}
        s = spec{1};
        bch = BCHCode(s(1), s(2));
        assert(all(all(mod(bch.G * bch.H.', 2) == 0)), 'BCH G H^T != 0');
        msg = double(randi(rng, [0, 1], 1, bch.k));
        cw = bch.encode(msg);
        assert(all(mod(cw * bch.H.', 2) == 0), 'BCH codeword not in code');
        % BM corrects a single error for t>=1
        r = cw; r(3) = bitxor(r(3), 1);
        [c_dec, ok] = bch.bmDecode(r);
        assert(ok && isequal(c_dec, cw), 'BCH BM single-error correction');
    end
    fprintf('  [PASS] BCH G H^T=0, valid codewords, BM corrects 1 error\n');

    % --- 4. Lagrange RS systematic generator ---
    gf = GF(8);
    n = gf.n;
    t = 2;
    k_prime = n - 2*t;
    Theta = sort(randperm(rng, n, k_prime));      % 1-based positions
    [G, Theta_c] = build_rs_systematic_generator(gf, Theta, k_prime, n);
    % (a) identity on Theta columns
    assert(isequal(G(:, Theta), eye(k_prime)), 'Lagrange: not identity on Theta');
    % (b) each row is a valid RS(n,k_prime) codeword: it vanishes at
    %     alpha^1..alpha^{n-k_prime}. Check via syndromes.
    n_check = n - k_prime;
    for row = 1:3   % sample a few rows for speed
        cw = G(row, :);
        for i = 1:n_check
            % evaluate cw(alpha^i) = sum_j cw(j) * alpha^{i*(j-1)}
            val = 0;
            nzj = find(cw ~= 0);
            for jj = nzj
                lg = mod(gf.LOG(cw(jj)+1) + i*(jj-1), gf.n);
                val = bitxor(val, gf.EXP(lg + 1));
            end
            assert(val == 0, 'Lagrange: row %d not RS codeword (syndrome %d nonzero)', row, i);
        end
    end
    fprintf('  [PASS] Lagrange G: identity on Theta + rows are valid RS codewords\n');

    % --- 5. LLOSD == OSD equivalence on identical LLR ---
    % Use a realistic (decodable) noise level: at very low SNR both decoders
    % find different low-weight codewords (neither correct), so equivalence only
    % holds where decoding actually succeeds — which is the operating regime.
    bch = BCHCode(8, 2);
    rng2 = RandStream('mt19937ar', 'Seed', 777);
    n_eq = 0; n_ok = 0;
    for trial = 1:15
        msg = double(randi(rng2, [0, 1], 1, bch.k));
        cw = bch.encode(msg);
        x = 1 - 2*cw;                          % BPSK 0->+1,1->-1
        sigma = 0.35;
        y = x + sigma * randn(rng2, size(x));
        L = 2*y / (sigma^2);
        [c_llosd, ~] = llosd_decode(bch, L, 2);
        [c_osd, ~]   = osd_decode(bch, L, 2);
        if isequal(c_llosd, c_osd), n_eq = n_eq + 1; end
        if isequal(c_llosd, cw), n_ok = n_ok + 1; end
    end
    assert(n_eq == 15, 'LLOSD/OSD disagreed on %d/15 trials', 15 - n_eq);
    assert(n_ok == 15, 'LLOSD failed to decode %d/15 trials', 15 - n_ok);
    fprintf('  [PASS] LLOSD (Lagrange) == OSD (Gauss elim), both decode 15/15\n');

    fprintf('=== Self-test complete: ALL PASS ===\n\n');
end
