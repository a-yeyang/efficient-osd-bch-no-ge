classdef HardCascade
%HARDCASCADE MATLAB R2022b reproduction of Work 03's hard RS+BCH cascade.
%
% This class ports the executable logic in hc_src/bch_t2.py,
% hc_src/cascade.py, hc_src/latency_model.py, and hc_src/simulate.py.
% It intentionally owns its GF, RS, BCH, PAM4, and simulation primitives so
% that Work 03 can run from its own MATLAB directory with Base MATLAB only.
%
% Conventions retained from the Python reference:
%   * GF elements are integers in 0 .. 2^m-1 and use polynomial bases.
%   * Polynomial coefficient vectors are low-degree first.
%   * Systematic RS/BCH words are [parity, message].
%   * RS symbols are serialized LSB-first.
%   * PAM4 labels are 00->-3, 01->-1, 11->+1, 10->+3.

methods(Static)
    %% Paths and small utilities
    function p = paths()
        %PATHS Return stable, repository-independent paths for Work 03.
        matlabRoot = fileparts(mfilename('fullpath'));
        codeRoot = fileparts(matlabRoot);
        workRoot = fileparts(codeRoot);
        p = struct( ...
            'matlab_root', matlabRoot, ...
            'code_root', codeRoot, ...
            'work_root', workRoot, ...
            'assets_root', fullfile(workRoot, 'assets'), ...
            'data_root', fullfile(workRoot, 'assets', 'data'), ...
            'figures_root', fullfile(workRoot, 'assets', 'figures'), ...
            'logs_root', fullfile(workRoot, 'assets', 'logs'));
    end

    function p = ensure_asset_dirs()
        %ENSURE_ASSET_DIRS Create only Work 03's own output directories.
        p = HardCascade.paths();
        dirs = {p.assets_root, p.data_root, p.figures_root, p.logs_root};
        for i = 1:numel(dirs)
            if ~exist(dirs{i}, 'dir')
                mkdir(dirs{i});
            end
        end
    end

    function write_json(filename, value)
        %WRITE_JSON Write a JSON artifact with a parent-directory check.
        parent = fileparts(filename);
        if ~isempty(parent) && ~exist(parent, 'dir')
            mkdir(parent);
        end
        fid = fopen(filename, 'w');
        assert(fid >= 0, 'HardCascade:WriteJson', ...
            'Cannot open %s for writing.', filename);
        cleanup = onCleanup(@() fclose(fid));
        fprintf(fid, '%s\n', jsonencode(value, 'PrettyPrint', true));
    end

    function save_figure_pair(fig, outputBase)
        %SAVE_FIGURE_PAIR Save deterministic PNG/PDF artifacts.
        parent = fileparts(outputBase);
        if ~exist(parent, 'dir')
            mkdir(parent);
        end
        drawnow;
        exportgraphics(fig, [outputBase '.png'], 'Resolution', 140);
        exportgraphics(fig, [outputBase '.pdf'], 'ContentType', 'vector');
    end

    function counters = op_counters()
        %OP_COUNTERS Python-compatible operation-counter shape.
        counters = struct('f2m', 0);
    end

    %% GF(2^m) arithmetic (same primitive-polynomial table as Python)
    function gf = gf_init(m)
        %GF_INIT Construct EXP/LOG tables for 1 <= m <= 8.
        primitive = [0 3 7 11 19 37 67 137 285];
        validateattributes(m, {'numeric'}, {'scalar', 'integer', '>=', 1, '<=', 8});
        n = 2^m - 1;
        expTable = zeros(1, 2*n + 2);
        logTable = -ones(1, 2^m);
        x = 1;
        for exponent = 0:n-1
            expTable(exponent + 1) = x;
            logTable(x + 1) = exponent;
            x = bitshift(x, 1);
            if bitand(x, 2^m) ~= 0
                x = bitxor(x, primitive(m + 1));
            end
        end
        for exponent = n:2*n+1
            expTable(exponent + 1) = expTable(exponent - n + 1);
        end
        gf = struct('m', m, 'n', n, 'prim', primitive(m + 1), ...
            'EXP', expTable, 'LOG', logTable);
    end

    function c = gf_add(a, b)
        %GF_ADD Addition in characteristic two.
        c = bitxor(a, b);
    end

    function c = gf_mul(gf, a, b)
        %GF_MUL Scalar multiplication through EXP/LOG tables.
        assert(isscalar(a) && isscalar(b), ...
            'HardCascade:GfMulScalar', 'gf_mul accepts scalar operands.');
        a = double(a);
        b = double(b);
        if a == 0 || b == 0
            c = 0;
            return;
        end
        c = gf.EXP(gf.LOG(a + 1) + gf.LOG(b + 1) + 1);
    end

    function c = gf_div(gf, a, b)
        %GF_DIV Scalar division through EXP/LOG tables.
        assert(isscalar(a) && isscalar(b), ...
            'HardCascade:GfDivScalar', 'gf_div accepts scalar operands.');
        a = double(a);
        b = double(b);
        assert(b ~= 0, 'HardCascade:GfDivideByZero', 'Division by zero in GF.');
        if a == 0
            c = 0;
            return;
        end
        c = gf.EXP(mod(gf.LOG(a + 1) - gf.LOG(b + 1), gf.n) + 1);
    end

    function c = gf_inv(gf, a)
        %GF_INV Multiplicative inverse.
        assert(a ~= 0, 'HardCascade:GfInverseZero', 'Zero has no GF inverse.');
        c = gf.EXP(mod(gf.n - gf.LOG(double(a) + 1), gf.n) + 1);
    end

    function c = gf_pow(gf, a, exponent)
        %GF_POW Integer exponentiation in GF(2^m).
        if a == 0
            c = double(exponent == 0);
            return;
        end
        c = gf.EXP(mod(gf.LOG(double(a) + 1) * exponent, gf.n) + 1);
    end

    function out = poly_mul(gf, a, b)
        %POLY_MUL Low-degree-first polynomial multiplication over GF.
        a = double(a(:).');
        b = double(b(:).');
        out = zeros(1, numel(a) + numel(b) - 1);
        for i = 1:numel(a)
            if a(i) == 0
                continue;
            end
            for j = 1:numel(b)
                if b(j) ~= 0
                    out(i + j - 1) = bitxor(out(i + j - 1), ...
                        HardCascade.gf_mul(gf, a(i), b(j)));
                end
            end
        end
    end

    function [quotient, remainder] = poly_divmod(gf, numerator, denominator)
        %POLY_DIVMOD Low-degree-first long division over GF(2^m).
        numerator = double(numerator(:).');
        denominator = double(denominator(:).');
        while numel(denominator) > 1 && denominator(end) == 0
            denominator(end) = [];
        end
        assert(any(denominator), 'HardCascade:ZeroDivisor', 'Zero polynomial divisor.');
        remainder = numerator;
        quotient = zeros(1, max(0, numel(remainder) - numel(denominator) + 1));
        inverseLead = HardCascade.gf_inv(gf, denominator(end));
        for i = numel(quotient):-1:1
            leadIndex = i + numel(denominator) - 1;
            if remainder(leadIndex) == 0
                continue;
            end
            coefficient = HardCascade.gf_mul(gf, remainder(leadIndex), inverseLead);
            quotient(i) = coefficient;
            for j = 1:numel(denominator)
                remainder(i + j - 1) = bitxor(remainder(i + j - 1), ...
                    HardCascade.gf_mul(gf, coefficient, denominator(j)));
            end
        end
        while ~isempty(remainder) && remainder(end) == 0
            remainder(end) = [];
        end
    end

    %% BCH(t=2): generator, systematic encoding, conventional/direct decode
    function bch = bch_t2_create(m, buildGenerator)
        %BCH_T2_CREATE Create primitive narrow-sense BCH(2^m-1, n-2m, t=2).
        %
        % buildGenerator is optional and defaults to true.  It supplies the
        % Python class's systematic k-by-n G matrix for inspection/tests;
        % encoding itself uses the equivalent polynomial-division path.
        if nargin < 2
            buildGenerator = true;
        end
        validateattributes(m, {'numeric'}, {'scalar', 'integer', '>=', 3, '<=', 8});
        gf = HardCascade.gf_init(m);
        seen = [];
        generator = 1;
        for rootExponent = [1 3]
            coset = HardCascade.cyclotomic_coset(rootExponent, gf.n);
            representative = min(coset);
            if ismember(representative, seen)
                continue;
            end
            seen(end + 1) = representative; %#ok<AGROW>
            minimal = 1;
            for exponent = coset
                minimal = HardCascade.poly_mul(gf, minimal, [gf.EXP(exponent + 1), 1]);
            end
            assert(all(minimal == 0 | minimal == 1), ...
                'HardCascade:NonBinaryMinimalPolynomial', ...
                'A BCH minimal polynomial must have binary coefficients.');
            generator = HardCascade.poly_mul(gf, generator, minimal);
            generator = mod(generator, 2);
        end
        n = gf.n;
        k = n - (numel(generator) - 1);
        assert(k == n - 2*m, 'HardCascade:UnexpectedBchDimension', ...
            'BCH(t=2) dimension must equal n - 2m.');
        bch = struct('m', m, 'gf', gf, 'n', n, 'k', k, 't', 2, ...
            'd_design', 5, 'g_poly', generator, ...
            'lut', HardCascade.build_lut_A(gf), 'G', []);
        if buildGenerator
            G = zeros(k, n);
            for i = 1:k
                unitMessage = zeros(1, k);
                unitMessage(i) = 1;
                G(i, :) = HardCascade.bch_t2_encode(bch, unitMessage);
            end
            bch.G = G;
        end
    end

    function coset = cyclotomic_coset(exponent, n)
        %CYCLOTOMIC_COSET Return the binary cyclotomic coset modulo n.
        coset = [];
        value = mod(exponent, n);
        while ~ismember(value, coset)
            coset(end + 1) = value; %#ok<AGROW>
            value = mod(2 * value, n);
        end
    end

    function codeword = bch_t2_encode(bch, message)
        %BCH_T2_ENCODE Systematic [parity, message] BCH encoding.
        message = double(message(:).');
        assert(numel(message) == bch.k, 'HardCascade:BchMessageLength', ...
            'Expected a BCH message of length %d.', bch.k);
        assert(all(message == 0 | message == 1), 'HardCascade:BchBinaryMessage', ...
            'BCH messages must be binary.');
        parityLength = bch.n - bch.k;
        dividend = [zeros(1, parityLength), message];
        [~, remainder] = HardCascade.poly_divmod(bch.gf, dividend, bch.g_poly);
        parity = [remainder, zeros(1, parityLength - numel(remainder))];
        codeword = [parity, message];
        assert(numel(codeword) == bch.n, 'HardCascade:BchEncodeLength', ...
            'Unexpected BCH codeword length.');
    end

    function message = bch_extract_message(bch, codeword)
        %BCH_EXTRACT_MESSAGE Extract a BCH systematic message suffix.
        codeword = double(codeword(:).');
        assert(numel(codeword) == bch.n, 'HardCascade:BchWordLength', ...
            'Expected a BCH word of length %d.', bch.n);
        message = codeword(bch.n - bch.k + 1:end);
    end

    function [s1, s3, counters] = bch_syndromes(bch, received, counters)
        %BCH_SYNDROMES Compute the odd BCH syndromes S1 and S3.
        if nargin < 3 || isempty(counters)
            counters = HardCascade.op_counters();
        end
        received = double(received(:).');
        assert(numel(received) == bch.n && all(received == 0 | received == 1), ...
            'HardCascade:BchReceivedWord', 'Expected a binary BCH word of length n.');
        positions = find(received == 1) - 1;
        s1 = 0;
        s3 = 0;
        if isempty(positions)
            return;
        end
        for position = positions
            s1 = bitxor(s1, bch.gf.EXP(mod(position, bch.gf.n) + 1));
            s3 = bitxor(s3, bch.gf.EXP(mod(3 * position, bch.gf.n) + 1));
        end
        % bch_t2.py accounts for a fully-parallel 2n syndrome network.
        counters.f2m = counters.f2m + 2 * bch.n;
    end

    function lut = build_lut_A(gf)
        %BUILD_LUT_A LUT for roots of A(Y)=Y^2+Y+k over GF(2^m).
        fieldSize = 2^gf.m;
        lut = struct('roots', zeros(fieldSize, 2), 'valid', false(fieldSize, 1));
        for x = 0:fieldSize-1
            k = bitxor(HardCascade.gf_mul(gf, x, x), x);
            if ~lut.valid(k + 1)
                lut.roots(k + 1, :) = [x, bitxor(x, 1)];
                lut.valid(k + 1) = true;
            end
        end
    end

    function [decoded, ok, counters] = bch_decode_conventional(bch, received, counters)
        %BCH_DECODE_CONVENTIONAL t=2 BM-equivalent locator + Chien decoder.
        %
        % For a binary t=2 BCH code, S1/S3 fully determine the simplified
        % Berlekamp--Massey result.  The two-error branch performs the same
        % Chien root search as the conventional Python implementation.
        if nargin < 3 || isempty(counters)
            counters = HardCascade.op_counters();
        end
        received = double(received(:).');
        [s1, s3, counters] = HardCascade.bch_syndromes(bch, received, counters);
        s1Cubed = HardCascade.gf_mul(bch.gf, HardCascade.gf_mul(bch.gf, s1, s1), s1);
        counters.f2m = counters.f2m + 2;

        if s1 == 0 && s3 == 0
            decoded = received;
            ok = true;
            return;
        end
        if s1 == 0 && s3 ~= 0
            decoded = received;
            ok = false;
            return;
        end
        if s1Cubed == s3
            position = mod(bch.gf.LOG(s1 + 1), bch.n);
            decoded = received;
            decoded(position + 1) = 1 - decoded(position + 1);
            ok = true;
            return;
        end

        constantTerm = HardCascade.gf_div(bch.gf, bitxor(s1Cubed, s3), s1);
        counters.f2m = counters.f2m + 2;
        errorPositions = [];
        for position = 0:bch.n-1
            x = bch.gf.EXP(position + 1);
            xSquared = bch.gf.EXP(mod(2 * position, bch.n) + 1);
            value = bitxor(bitxor(xSquared, HardCascade.gf_mul(bch.gf, s1, x)), constantTerm);
            counters.f2m = counters.f2m + 2;
            if value == 0
                errorPositions(end + 1) = position; %#ok<AGROW>
                if numel(errorPositions) == 2
                    break;
                end
            end
        end
        if numel(errorPositions) ~= 2
            decoded = received;
            ok = false;
            return;
        end
        decoded = received;
        decoded(errorPositions + 1) = 1 - decoded(errorPositions + 1);
        ok = true;
    end

    function [decoded, ok, counters] = bch_decode_direct(bch, received, counters)
        %BCH_DECODE_DIRECT LUT direct-root decoder from bch_t2.py.
        if nargin < 3 || isempty(counters)
            counters = HardCascade.op_counters();
        end
        received = double(received(:).');
        [s1, s3, counters] = HardCascade.bch_syndromes(bch, received, counters);
        s1Cubed = HardCascade.gf_mul(bch.gf, HardCascade.gf_mul(bch.gf, s1, s1), s1);
        counters.f2m = counters.f2m + 2;

        if s1 == 0 && s3 == 0
            decoded = received;
            ok = true;
            return;
        end
        if s1 == 0 && s3 ~= 0
            decoded = received;
            ok = false;
            return;
        end
        if s1Cubed == s3
            position = mod(bch.gf.LOG(s1 + 1), bch.n);
            decoded = received;
            decoded(position + 1) = 1 - decoded(position + 1);
            ok = true;
            return;
        end

        lutIndex = HardCascade.gf_div(bch.gf, bitxor(s1Cubed, s3), s1Cubed);
        counters.f2m = counters.f2m + 1;
        if ~bch.lut.valid(lutIndex + 1)
            decoded = received;
            ok = false;
            return;
        end
        y1 = bch.lut.roots(lutIndex + 1, 1);
        y2 = bch.lut.roots(lutIndex + 1, 2);
        x1 = HardCascade.gf_mul(bch.gf, s1, y1);
        x2 = HardCascade.gf_mul(bch.gf, s1, y2);
        counters.f2m = counters.f2m + 2;
        if x1 == 0 || x2 == 0
            decoded = received;
            ok = false;
            return;
        end
        position1 = mod(bch.gf.LOG(x1 + 1), bch.n);
        position2 = mod(bch.gf.LOG(x2 + 1), bch.n);
        if position1 == position2
            decoded = received;
            ok = false;
            return;
        end
        decoded = received;
        decoded([position1, position2] + 1) = 1 - decoded([position1, position2] + 1);
        ok = true;
    end

    %% Reed--Solomon code (hard BM + Chien + Forney)
    function rs = rs_create(m, k)
        %RS_CREATE Primitive narrow-sense systematic RS(n,k) code.
        gf = HardCascade.gf_init(m);
        n = gf.n;
        validateattributes(k, {'numeric'}, {'scalar', 'integer', '>=', 1, '<=', n});
        assert(mod(n - k, 2) == 0, 'HardCascade:RsParity', ...
            'This BM construction requires an even RS parity-symbol count.');
        t = (n - k) / 2;
        generator = 1;
        for i = 1:2*t
            generator = HardCascade.poly_mul(gf, generator, [gf.EXP(i + 1), 1]);
        end
        assert(numel(generator) - 1 == n - k, 'HardCascade:RsGeneratorDegree', ...
            'Unexpected RS generator degree.');
        rs = struct('gf', gf, 'm', m, 'n', n, 'k', k, 't', t, ...
            'd', n - k + 1, 'g_poly', generator, ...
            'alpha_pow', gf.EXP(1:n));
    end

    function codeword = rs_encode_systematic(rs, message)
        %RS_ENCODE_SYSTEMATIC Systematic [parity, message] RS encoding.
        message = double(message(:).');
        assert(numel(message) == rs.k, 'HardCascade:RsMessageLength', ...
            'Expected an RS message of length %d.', rs.k);
        assert(all(message >= 0 & message < 2^rs.m & message == floor(message)), ...
            'HardCascade:RsMessageAlphabet', 'RS symbols must belong to GF(2^m).');
        parityLength = rs.n - rs.k;
        [~, remainder] = HardCascade.poly_divmod(rs.gf, ...
            [zeros(1, parityLength), message], rs.g_poly);
        parity = [remainder, zeros(1, parityLength - numel(remainder))];
        codeword = [parity, message];
    end

    function message = rs_extract_message(rs, codeword)
        %RS_EXTRACT_MESSAGE Extract systematic RS message symbols.
        codeword = double(codeword(:).');
        assert(numel(codeword) == rs.n, 'HardCascade:RsWordLength', ...
            'Expected an RS word of length %d.', rs.n);
        message = codeword(rs.n - rs.k + 1:end);
    end

    function [decoded, ok, counters] = rs_bm_decode(rs, received, counters)
        %RS_BM_DECODE Port of RSCode.bm_decode from Work 02's Python source.
        if nargin < 3 || isempty(counters)
            counters = HardCascade.op_counters();
        end
        received = double(received(:).');
        assert(numel(received) == rs.n, 'HardCascade:RsReceivedLength', ...
            'Expected an RS word of length %d.', rs.n);
        assert(all(received >= 0 & received < 2^rs.m & received == floor(received)), ...
            'HardCascade:RsReceivedAlphabet', 'Received RS symbols must belong to GF(2^m).');

        gf = rs.gf;
        nonzeroPositions = find(received ~= 0) - 1;
        syndromes = zeros(1, 2 * rs.t + 1);
        if isempty(nonzeroPositions)
            decoded = received;
            ok = true;
            return;
        end
        for syndromeIndex = 1:2*rs.t
            syndrome = 0;
            for position = nonzeroPositions
                term = HardCascade.gf_mul(gf, received(position + 1), ...
                    gf.EXP(mod(syndromeIndex * position, gf.n) + 1));
                syndrome = bitxor(syndrome, term);
            end
            syndromes(syndromeIndex + 1) = syndrome;
        end
        counters.f2m = counters.f2m + numel(nonzeroPositions) * 2 * rs.t;
        if ~any(syndromes(2:end))
            decoded = received;
            ok = true;
            return;
        end

        % Berlekamp--Massey recurrence.
        locator = 1;
        auxiliary = 1;
        discrepancyScale = 1;
        degree = 0;
        shift = 1;
        for iteration = 1:2*rs.t
            discrepancy = syndromes(iteration + 1);
            for i = 1:degree
                if i + 1 <= numel(locator) && locator(i + 1) ~= 0
                    discrepancy = bitxor(discrepancy, ...
                        HardCascade.gf_mul(gf, locator(i + 1), syndromes(iteration - i + 1)));
                    counters.f2m = counters.f2m + 1;
                end
            end
            if discrepancy == 0
                shift = shift + 1;
                continue;
            end
            coefficient = HardCascade.gf_div(gf, discrepancy, discrepancyScale);
            counters.f2m = counters.f2m + 1;
            shiftedAuxiliary = [zeros(1, shift), auxiliary];
            newLength = max(numel(locator), numel(shiftedAuxiliary));
            candidate = [locator, zeros(1, newLength - numel(locator))];
            paddedAuxiliary = zeros(1, newLength);
            paddedAuxiliary(1:numel(shiftedAuxiliary)) = shiftedAuxiliary;
            shiftedAuxiliary = paddedAuxiliary;
            for i = 1:newLength
                candidate(i) = bitxor(candidate(i), ...
                    HardCascade.gf_mul(gf, coefficient, shiftedAuxiliary(i)));
                counters.f2m = counters.f2m + 1;
            end
            if 2 * degree <= iteration - 1
                newDegree = iteration - degree;
                auxiliary = locator;
                discrepancyScale = discrepancy;
                locator = candidate;
                degree = newDegree;
                shift = 1;
            else
                locator = candidate;
                shift = shift + 1;
            end
        end

        % Chien search: roots Lambda(alpha^{-p}) identify symbol positions p.
        nonzeroDegrees = find(locator ~= 0) - 1;
        if isempty(nonzeroDegrees)
            decoded = received;
            ok = false;
            return;
        end
        errorPositions = [];
        for position = 0:rs.n-1
            value = 0;
            for degreeIndex = nonzeroDegrees
                value = bitxor(value, HardCascade.gf_mul(gf, locator(degreeIndex + 1), ...
                    gf.EXP(mod((rs.n - position) * degreeIndex, gf.n) + 1)));
            end
            if value == 0
                errorPositions(end + 1) = position; %#ok<AGROW>
            end
        end
        counters.f2m = counters.f2m + rs.n * numel(nonzeroDegrees);
        if numel(errorPositions) ~= degree || degree > rs.t
            decoded = received;
            ok = false;
            return;
        end

        % Forney error evaluation, matching the reference's syndrome origin.
        syndromePolynomial = [0, syndromes(2:end)];
        omega = zeros(1, 2*rs.t + 1);
        for i = 0:numel(syndromePolynomial)-1
            if syndromePolynomial(i + 1) == 0
                continue;
            end
            for j = 0:numel(locator)-1
                if i + j <= 2*rs.t && locator(j + 1) ~= 0
                    omega(i + j + 1) = bitxor(omega(i + j + 1), ...
                        HardCascade.gf_mul(gf, syndromePolynomial(i + 1), locator(j + 1)));
                    counters.f2m = counters.f2m + 1;
                end
            end
        end
        locatorDerivative = zeros(1, numel(locator));
        for degreeIndex = 1:numel(locator)-1
            if mod(degreeIndex, 2) == 1
                locatorDerivative(degreeIndex) = locator(degreeIndex + 1);
            end
        end

        decoded = received;
        for position = errorPositions
            omegaValue = 0;
            derivativeValue = 0;
            for i = 0:numel(omega)-1
                if omega(i + 1) ~= 0
                    omegaValue = bitxor(omegaValue, HardCascade.gf_mul(gf, omega(i + 1), ...
                        gf.EXP(mod(i * (rs.n - position), gf.n) + 1)));
                    counters.f2m = counters.f2m + 1;
                end
            end
            for i = 0:numel(locatorDerivative)-1
                if locatorDerivative(i + 1) ~= 0
                    derivativeValue = bitxor(derivativeValue, ...
                        HardCascade.gf_mul(gf, locatorDerivative(i + 1), ...
                        gf.EXP(mod(i * (rs.n - position), gf.n) + 1)));
                    counters.f2m = counters.f2m + 1;
                end
            end
            if derivativeValue == 0
                decoded = received;
                ok = false;
                return;
            end
            errorValue = HardCascade.gf_div(gf, omegaValue, derivativeValue);
            errorValue = HardCascade.gf_mul(gf, gf.EXP(position + 1), errorValue);
            counters.f2m = counters.f2m + 1;
            decoded(position + 1) = bitxor(decoded(position + 1), errorValue);
        end
        ok = true;
    end

    %% Serialization and PAM4 hard-decision channel
    function bits = symbols_to_bits(symbols, m)
        %SYMBOLS_TO_BITS Serialize GF symbols LSB-first, like Python.
        symbols = double(symbols(:).');
        assert(all(symbols >= 0 & symbols < 2^m & symbols == floor(symbols)), ...
            'HardCascade:SymbolAlphabet', 'Symbols must belong to GF(2^m).');
        bits = zeros(1, numel(symbols) * m);
        for i = 1:numel(symbols)
            for bit = 1:m
                bits((i - 1) * m + bit) = bitget(symbols(i), bit);
            end
        end
    end

    function symbols = bits_to_symbols(bits, m)
        %BITS_TO_SYMBOLS Parse LSB-first bits into GF symbols.
        bits = double(bits(:).');
        assert(mod(numel(bits), m) == 0 && all(bits == 0 | bits == 1), ...
            'HardCascade:SymbolBits', 'Expected a binary stream divisible by m.');
        symbols = zeros(1, numel(bits) / m);
        for i = 1:numel(symbols)
            for bit = 1:m
                symbols(i) = symbols(i) + bits((i - 1) * m + bit) * 2^(bit - 1);
            end
        end
    end

    function channelSymbols = bits_to_pam4(bits)
        %BITS_TO_PAM4 Gray map pairs [b1,b0] to {-3,-1,+1,+3}.
        bits = double(bits(:).');
        assert(mod(numel(bits), 2) == 0 && all(bits == 0 | bits == 1), ...
            'HardCascade:Pam4Bits', 'PAM4 input must be an even-length binary vector.');
        channelSymbols = zeros(1, numel(bits) / 2);
        for i = 1:numel(channelSymbols)
            b1 = bits(2*i - 1);
            b0 = bits(2*i);
            if b1 == 0 && b0 == 0
                channelSymbols(i) = -3;
            elseif b1 == 0 && b0 == 1
                channelSymbols(i) = -1;
            elseif b1 == 1 && b0 == 1
                channelSymbols(i) = 1;
            else
                channelSymbols(i) = 3;
            end
        end
    end

    function bits = pam4_to_bits_hard(received)
        %PAM4_TO_BITS_HARD Nearest-neighbour PAM4 hard demodulation.
        received = double(received(:).');
        levels = [-3, -1, 1, 3];
        labels = [0 0; 0 1; 1 1; 1 0];
        bits = zeros(1, 2 * numel(received));
        for i = 1:numel(received)
            [~, index] = min(abs(levels - received(i)));
            bits(2*i - 1:2*i) = labels(index, :);
        end
    end

    function sigma = sigma_from_ebn0_pam4(ebn0Db, rate)
        %SIGMA_FROM_EBN0_PAM4 Exact Work 03/Python PAM4 AWGN convention.
        % sigma^2 = E_s * rate / (4 * 10^(Eb/N0/10)), E_s = 5.
        validateattributes(rate, {'numeric'}, {'scalar', 'positive'});
        sigma = sqrt(5 * rate / (4 * 10^(double(ebn0Db) / 10)));
    end

    function received = awgn_channel(transmitted, sigma, source)
        %AWGN_CHANNEL Add Gaussian noise; source may be a RandStream or seed.
        transmitted = double(transmitted(:).');
        if nargin < 3 || isempty(source)
            noise = randn(size(transmitted));
        elseif isa(source, 'RandStream')
            noise = randn(source, size(transmitted));
        elseif isnumeric(source) && isscalar(source)
            stream = RandStream('mt19937ar', 'Seed', double(source));
            noise = randn(stream, size(transmitted));
        else
            error('HardCascade:NoiseSource', ...
                'source must be empty, a RandStream, or a numeric scalar seed.');
        end
        received = transmitted + sigma * noise;
    end

    function hardBits = run_pam4_channel_hard(bits, ebn0Db, rate, source)
        %RUN_PAM4_CHANNEL_HARD PAM4 modulation, AWGN, and hard demodulation.
        if nargin < 4
            source = [];
        end
        bits = double(bits(:).');
        assert(all(bits == 0 | bits == 1), 'HardCascade:ChannelBits', ...
            'Channel input must be binary.');
        originalLength = numel(bits);
        if mod(originalLength, 2) ~= 0
            mappedBits = [bits, 0];
        else
            mappedBits = bits;
        end
        sigma = HardCascade.sigma_from_ebn0_pam4(ebn0Db, rate);
        received = HardCascade.awgn_channel(HardCascade.bits_to_pam4(mappedBits), sigma, source);
        hardBits = HardCascade.pam4_to_bits_hard(received);
        hardBits = hardBits(1:originalLength);
    end

    %% Hard cascade and pure-RS codecs
    function cfg = hard_cascade_config(m, kRs)
        %HARD_CASCADE_CONFIG MATLAB counterpart of HardCascadeConfig.
        nRs = 2^m - 1;
        validateattributes(kRs, {'numeric'}, {'scalar', 'integer', '>=', 1, '<=', nRs});
        assert(mod(nRs - kRs, 2) == 0, 'HardCascade:RsConfigParity', ...
            'k_rs must give an integer RS correction radius.');
        cfg = struct('m', m, 'k_rs', kRs, 'n_rs', nRs, ...
            't_rs', (nRs - kRs) / 2);
        cfg.description = sprintf('RS(%d,%d, t=%d) + BCH(%d, t=2)', ...
            cfg.n_rs, cfg.k_rs, cfg.t_rs, cfg.n_rs);
    end

    function codec = cascade_create(m, kRs, decoder)
        %CASCADE_CREATE Create a hard RS+BCH cascade ('conv' or 'direct').
        if nargin < 3 || isempty(decoder)
            decoder = 'direct';
        end
        decoder = lower(char(decoder));
        assert(ismember(decoder, {'conv', 'direct'}), 'HardCascade:Decoder', ...
            'decoder must be ''conv'' or ''direct''.');
        cfg = HardCascade.hard_cascade_config(m, kRs);
        rs = HardCascade.rs_create(m, kRs);
        bch = HardCascade.bch_t2_create(m, true);
        rsBits = m * rs.n;
        paddingBits = mod(-rsBits, bch.k);
        blockCount = (rsBits + paddingBits) / bch.k;
        codec = struct('kind', 'cascade', 'cfg', cfg, 'm', m, 'k_rs', kRs, ...
            'rs', rs, 'bch', bch, 'decoder', decoder, ...
            'rs_bits', rsBits, 'n_pad_bits', paddingBits, ...
            'n_bch_blocks', blockCount, 'n_coded_bits', blockCount * bch.n, ...
            'n_info_bits', kRs * m, ...
            'effective_rate', (kRs * m) / (blockCount * bch.n));
    end

    function codec = pure_rs_create(m, kRs)
        %PURE_RS_CREATE Pure hard-decision RS baseline codec.
        cfg = HardCascade.hard_cascade_config(m, kRs);
        rs = HardCascade.rs_create(m, kRs);
        codedBits = m * rs.n;
        codec = struct('kind', 'pure_rs', 'cfg', cfg, 'm', m, 'k_rs', kRs, ...
            'rs', rs, 'rs_bits', codedBits, 'n_coded_bits', codedBits, ...
            'n_info_bits', kRs * m, 'effective_rate', (kRs * m) / codedBits);
    end

    function bits = cascade_encode(codec, messageSymbols)
        %CASCADE_ENCODE RS encode -> LSB-first bits -> BCH blocks.
        assert(isfield(codec, 'kind') && strcmp(codec.kind, 'cascade'), ...
            'HardCascade:CodecKind', 'Expected a cascade codec.');
        messageSymbols = double(messageSymbols(:).');
        assert(numel(messageSymbols) == codec.k_rs, ...
            'HardCascade:CascadeMessageLength', 'Wrong cascade message length.');
        rsCodeword = HardCascade.rs_encode_systematic(codec.rs, messageSymbols);
        rawBits = HardCascade.symbols_to_bits(rsCodeword, codec.m);
        if codec.n_pad_bits > 0
            rawBits = [rawBits, zeros(1, codec.n_pad_bits)];
        end
        bits = zeros(1, codec.n_coded_bits);
        for block = 1:codec.n_bch_blocks
            inputIndices = (block - 1) * codec.bch.k + 1:block * codec.bch.k;
            outputIndices = (block - 1) * codec.bch.n + 1:block * codec.bch.n;
            bits(outputIndices) = HardCascade.bch_t2_encode(codec.bch, rawBits(inputIndices));
        end
    end

    function [messageSymbols, okRs, stats] = cascade_decode(codec, hardBits, counters)
        %CASCADE_DECODE Hard inner BCH decoding followed by RS BM decoding.
        if nargin < 3 || isempty(counters)
            counters = HardCascade.op_counters();
        end
        assert(isfield(codec, 'kind') && strcmp(codec.kind, 'cascade'), ...
            'HardCascade:CodecKind', 'Expected a cascade codec.');
        hardBits = double(hardBits(:).');
        assert(numel(hardBits) == codec.n_coded_bits && all(hardBits == 0 | hardBits == 1), ...
            'HardCascade:CascadeHardBits', 'Wrong cascade hard-bit vector.');
        recoveredBits = zeros(1, codec.n_bch_blocks * codec.bch.k);
        bchSuccesses = 0;
        for block = 1:codec.n_bch_blocks
            wordIndices = (block - 1) * codec.bch.n + 1:block * codec.bch.n;
            if strcmp(codec.decoder, 'conv')
                [decodedWord, bchOk, counters] = ...
                    HardCascade.bch_decode_conventional(codec.bch, hardBits(wordIndices), counters);
            else
                [decodedWord, bchOk, counters] = ...
                    HardCascade.bch_decode_direct(codec.bch, hardBits(wordIndices), counters);
            end
            if bchOk
                bchSuccesses = bchSuccesses + 1;
            end
            messageIndices = (block - 1) * codec.bch.k + 1:block * codec.bch.k;
            recoveredBits(messageIndices) = HardCascade.bch_extract_message(codec.bch, decodedWord);
        end
        rsReceived = HardCascade.bits_to_symbols(recoveredBits(1:codec.rs_bits), codec.m);
        [decodedRs, okRs, counters] = HardCascade.rs_bm_decode(codec.rs, rsReceived, counters);
        if okRs
            messageSymbols = HardCascade.rs_extract_message(codec.rs, decodedRs);
        else
            messageSymbols = rsReceived(codec.rs.n - codec.rs.k + 1:end);
        end
        stats = struct('n_bch_ok', bchSuccesses, 'counters', counters);
    end

    function bits = pure_rs_encode(codec, messageSymbols)
        %PURE_RS_ENCODE RS encode and serialize LSB-first.
        assert(isfield(codec, 'kind') && strcmp(codec.kind, 'pure_rs'), ...
            'HardCascade:CodecKind', 'Expected a pure-RS codec.');
        bits = HardCascade.symbols_to_bits( ...
            HardCascade.rs_encode_systematic(codec.rs, messageSymbols), codec.m);
    end

    function [messageSymbols, ok, stats] = pure_rs_decode(codec, hardBits, counters)
        %PURE_RS_DECODE Hard-bit deserialization followed by RS BM decoding.
        if nargin < 3 || isempty(counters)
            counters = HardCascade.op_counters();
        end
        assert(isfield(codec, 'kind') && strcmp(codec.kind, 'pure_rs'), ...
            'HardCascade:CodecKind', 'Expected a pure-RS codec.');
        hardBits = double(hardBits(:).');
        assert(numel(hardBits) == codec.n_coded_bits && all(hardBits == 0 | hardBits == 1), ...
            'HardCascade:PureRsHardBits', 'Wrong pure-RS hard-bit vector.');
        received = HardCascade.bits_to_symbols(hardBits, codec.m);
        [decoded, ok, counters] = HardCascade.rs_bm_decode(codec.rs, received, counters);
        if ok
            messageSymbols = HardCascade.rs_extract_message(codec.rs, decoded);
        else
            messageSymbols = received(codec.rs.n - codec.rs.k + 1:end);
        end
        stats = struct('counters', counters);
    end

    %% Lagendijk-style latency model
    function cycles = rs_bm_cycles(t)
        %RS_BM_CYCLES 1 syndrome + 2t BM + Chien + Forney.
        cycles = 2 * t + 3;
    end

    function cycles = bch_conv_cycles(t, m) %#ok<INUSD>
        %BCH_CONV_CYCLES Table-VI/extrapolated conventional BCH latency.
        if t == 2
            cycles = 8;
        elseif t == 3
            cycles = 8;
        elseif t == 4
            cycles = 10;
        else
            cycles = 2 * t + 4;
        end
    end

    function cycles = bch_direct_cycles(t, m)
        %BCH_DIRECT_CYCLES Direct-root latency, including GF(2^7) refinement.
        if nargin < 2
            m = 8;
        end
        if t == 2
            if m <= 7
                cycles = 2;
            else
                cycles = 3;
            end
        elseif t == 3
            cycles = 4;
        elseif t == 4
            cycles = 8;
        else
            cycles = NaN;
        end
    end

    function cycles = cascade_serial(bchCycles, rsT)
        cycles = bchCycles + HardCascade.rs_bm_cycles(rsT);
    end

    function cycles = cascade_lagrange_v1(bchCycles, rsT)
        cycles = bchCycles + HardCascade.rs_bm_cycles(rsT) - 1;
    end

    function cycles = cascade_lagrange_v2(bchCycles, rsT)
        cycles = bchCycles + HardCascade.rs_bm_cycles(rsT) - 2;
    end

    function cycles = cascade_lagrange_shared(bchCycles, rsT)
        %CASCADE_LAGRANGE_SHARED Backward-compatible v1 alias.
        cycles = HardCascade.cascade_lagrange_v1(bchCycles, rsT);
    end

    function cycles = latency_cycles(codec, mode)
        %LATENCY_CYCLES Total cascade latency for none/v1/v2 sharing.
        if nargin < 2 || isempty(mode)
            mode = 'none';
        end
        assert(isfield(codec, 'kind') && strcmp(codec.kind, 'cascade'), ...
            'HardCascade:CodecKind', 'latency_cycles requires a cascade codec.');
        if strcmp(codec.decoder, 'conv')
            bchCycles = HardCascade.bch_conv_cycles(2, codec.m);
        else
            bchCycles = HardCascade.bch_direct_cycles(2, codec.m);
        end
        if islogical(mode)
            if mode
                mode = 'v1';
            else
                mode = 'none';
            end
        else
            mode = lower(char(mode));
        end
        switch mode
            case 'none'
                cycles = HardCascade.cascade_serial(bchCycles, codec.rs.t);
            case 'v1'
                cycles = HardCascade.cascade_lagrange_v1(bchCycles, codec.rs.t);
            case 'v2'
                cycles = HardCascade.cascade_lagrange_v2(bchCycles, codec.rs.t);
            otherwise
                error('HardCascade:LatencyMode', 'Unknown Lagrange mode: %s', mode);
        end
    end

    function [ratioPct, passes] = kpi_check(cascadeCycles, baselineCycles)
        %KPI_CHECK Return the +10%% latency KPI decision.
        ratioPct = (cascadeCycles - baselineCycles) / baselineCycles * 100;
        passes = ratioPct <= 10.0;
    end

    function summary = latency_summary_v1(cfg)
        %LATENCY_SUMMARY_V1 Exact main.py latency behavior.
        %
        % main.py intentionally calls bch_direct_cycles(2) without cfg.m;
        % preserving that default reproduces its n=127 v1 value of 20 cycles.
        rsOnly = HardCascade.rs_bm_cycles(cfg.t_rs);
        conventional = HardCascade.bch_conv_cycles(2);
        direct = HardCascade.bch_direct_cycles(2);
        summary = struct('pure_rs', rsOnly, ...
            'conv_no_share', HardCascade.cascade_serial(conventional, cfg.t_rs), ...
            'direct_no_share', HardCascade.cascade_serial(direct, cfg.t_rs), ...
            'conv_lagrange', HardCascade.cascade_lagrange_shared(conventional, cfg.t_rs), ...
            'direct_lagrange', HardCascade.cascade_lagrange_shared(direct, cfg.t_rs), ...
            'kpi_target', floor(rsOnly * 1.10));
    end

    function summary = latency_summary_v2(cfg)
        %LATENCY_SUMMARY_V2 Exact main_v2.py latency behavior.
        rsOnly = HardCascade.rs_bm_cycles(cfg.t_rs);
        conventional = HardCascade.bch_conv_cycles(2, cfg.m);
        direct = HardCascade.bch_direct_cycles(2, cfg.m);
        summary = struct('pure_rs', rsOnly, 'kpi_ceiling', rsOnly * 1.10, ...
            'bch_conv_cyc', conventional, 'bch_direct_cyc', direct, ...
            'conv_none', HardCascade.cascade_serial(conventional, cfg.t_rs), ...
            'conv_v1', HardCascade.cascade_lagrange_v1(conventional, cfg.t_rs), ...
            'conv_v2', HardCascade.cascade_lagrange_v2(conventional, cfg.t_rs), ...
            'direct_none', HardCascade.cascade_serial(direct, cfg.t_rs), ...
            'direct_v1', HardCascade.cascade_lagrange_v1(direct, cfg.t_rs), ...
            'direct_v2', HardCascade.cascade_lagrange_v2(direct, cfg.t_rs));
    end

    %% Monte-Carlo simulation driver
    function result = run_bench(methodName, codec, ebn0List, nInfoSymbols, mBits, ...
            seed, minFrameErrors, maxFrames, verbose)
        %RUN_BENCH MATLAB equivalent of hc_src.simulate.run_bench.
        if nargin < 6 || isempty(seed), seed = 0; end
        if nargin < 7 || isempty(minFrameErrors), minFrameErrors = 30; end
        if nargin < 8 || isempty(maxFrames), maxFrames = 2000; end
        if nargin < 9 || isempty(verbose), verbose = true; end
        ebn0List = double(ebn0List(:).');
        result = struct('ebn0_db', [], 'ser', [], 'ber', [], 'fer', [], ...
            'n_frames', [], 'n_frame_errors', [], 'avg_f2m_ops', []);
        for point = 1:numel(ebn0List)
            ebn0 = ebn0List(point);
            stream = RandStream('mt19937ar', 'Seed', seed + round(ebn0 * 100));
            frameCount = 0;
            frameErrors = 0;
            bitErrors = 0;
            symbolErrors = 0;
            accumulatedF2m = 0;
            started = tic;
            while frameCount < maxFrames
                message = randi(stream, 2^mBits, 1, nInfoSymbols) - 1;
                if strcmp(codec.kind, 'cascade')
                    coded = HardCascade.cascade_encode(codec, message);
                else
                    coded = HardCascade.pure_rs_encode(codec, message);
                end
                hard = HardCascade.run_pam4_channel_hard(coded, ebn0, codec.effective_rate, stream);
                counters = HardCascade.op_counters();
                if strcmp(codec.kind, 'cascade')
                    [messageHat, ~, stats] = HardCascade.cascade_decode(codec, hard, counters);
                else
                    [messageHat, ~, stats] = HardCascade.pure_rs_decode(codec, hard, counters);
                end
                frameCount = frameCount + 1;
                thisSymbolErrors = sum(messageHat ~= message);
                if thisSymbolErrors > 0
                    frameErrors = frameErrors + 1;
                end
                symbolErrors = symbolErrors + thisSymbolErrors;
                bitErrors = bitErrors + sum(HardCascade.symbols_to_bits(message, mBits) ~= ...
                    HardCascade.symbols_to_bits(messageHat, mBits));
                accumulatedF2m = accumulatedF2m + stats.counters.f2m;
                if frameErrors >= minFrameErrors && frameCount >= 100
                    break;
                end
            end
            ser = symbolErrors / max(1, frameCount * nInfoSymbols);
            ber = bitErrors / max(1, frameCount * nInfoSymbols * mBits);
            fer = frameErrors / max(1, frameCount);
            avgF2m = accumulatedF2m / max(1, frameCount);
            if verbose
                fprintf('  %s @ %.2f dB: FER=%.3e, BER=%.3e, avg_f2m=%.0f, %d frames, %.1fs\n', ...
                    methodName, ebn0, fer, ber, avgF2m, frameCount, toc(started));
            end
            result.ebn0_db(end + 1) = ebn0;
            result.ser(end + 1) = ser;
            result.ber(end + 1) = ber;
            result.fer(end + 1) = fer;
            result.n_frames(end + 1) = frameCount;
            result.n_frame_errors(end + 1) = frameErrors;
            result.avg_f2m_ops(end + 1) = avgF2m;
            if fer < 1e-6 && frameErrors < 2
                break;
            end
        end
    end

    %% Deterministic self-test used by the repository-level acceptance runner
    function report = selftest()
        %SELFTEST Exercise every numerical subsystem with real assertions.
        rng(20260723, 'twister');
        bch = HardCascade.bch_t2_create(5, true);
        nonzeroMessage = mod(0:bch.k-1, 2);
        codeword = HardCascade.bch_t2_encode(bch, nonzeroMessage);
        assert(isequal(codeword, mod(nonzeroMessage * bch.G, 2)), ...
            'HardCascade:GeneratorMatrix', 'Systematic BCH G does not encode identically.');
        [s1, s3] = HardCascade.bch_syndromes(bch, codeword);
        assert(s1 == 0 && s3 == 0 && any(codeword), ...
            'HardCascade:BchCodeword', 'Nonzero encoded BCH word is invalid.');

        decoderNames = {'conventional', 'direct'};
        decoderPass = false(1, 2);
        for decoderIndex = 1:2
            for errors = 0:2
                received = codeword;
                if errors > 0
                    positions = randperm(bch.n, errors);
                    received(positions) = 1 - received(positions);
                end
                if decoderIndex == 1
                    [decoded, ok] = HardCascade.bch_decode_conventional(bch, received);
                else
                    [decoded, ok] = HardCascade.bch_decode_direct(bch, received);
                end
                assert(ok && isequal(decoded, codeword), 'HardCascade:BchBoundedDistance', ...
                    '%s BCH decoder failed a %d-error word.', decoderNames{decoderIndex}, errors);
            end
            decoderPass(decoderIndex) = true;
        end

        % A three-error pattern with locators 1, alpha, 1+alpha has S1=0,
        % S3~=0 and is explicitly classified as uncorrectable by Table I.
        alpha = bch.gf.EXP(2);
        thirdLocator = bitxor(1, alpha);
        triple = [0, 1, bch.gf.LOG(thirdLocator + 1)];
        received = codeword;
        received(triple + 1) = 1 - received(triple + 1);
        [~, conventionalOk] = HardCascade.bch_decode_conventional(bch, received);
        [~, directOk] = HardCascade.bch_decode_direct(bch, received);
        assert(~conventionalOk && ~directOk, 'HardCascade:BchUncorrectable', ...
            'Known S1=0/S3~=0 pattern must be reported uncorrectable.');

        rs = HardCascade.rs_create(5, 27);
        rsMessage = mod(0:rs.k-1, 2^rs.m);
        rsWord = HardCascade.rs_encode_systematic(rs, rsMessage);
        [rsDecoded, rsOk] = HardCascade.rs_bm_decode(rs, rsWord);
        assert(rsOk && isequal(rsDecoded, rsWord), 'HardCascade:RsNoiseless', ...
            'Noiseless RS decoding failed.');
        rsReceived = rsWord;
        rsReceived([2, 11]) = bitxor(rsReceived([2, 11]), [3, 7]);
        [rsDecoded, rsOk] = HardCascade.rs_bm_decode(rs, rsReceived);
        assert(rsOk && isequal(rsDecoded, rsWord), 'HardCascade:RsBoundedDistance', ...
            'RS BM did not correct two symbol errors.');

        directCodec = HardCascade.cascade_create(5, 27, 'direct');
        cascadeBits = HardCascade.cascade_encode(directCodec, rsMessage);
        damagedBits = cascadeBits;
        damagedBits([2, 7]) = 1 - damagedBits([2, 7]);
        [decodedMessage, cascadeOk, stats] = HardCascade.cascade_decode(directCodec, damagedBits);
        assert(cascadeOk && isequal(decodedMessage, rsMessage) && stats.n_bch_ok == directCodec.n_bch_blocks, ...
            'HardCascade:CascadeBoundedDistance', 'Direct hard cascade failed a two-bit inner error.');
        conventionalCodec = HardCascade.cascade_create(5, 27, 'conv');
        [decodedMessage, cascadeOk] = HardCascade.cascade_decode(conventionalCodec, ...
            HardCascade.cascade_encode(conventionalCodec, rsMessage));
        assert(cascadeOk && isequal(decodedMessage, rsMessage), ...
            'HardCascade:CascadeConventional', 'Conventional hard cascade failed noiseless decoding.');

        pureCodec = HardCascade.pure_rs_create(5, 27);
        pureBits = HardCascade.pure_rs_encode(pureCodec, rsMessage);
        [decodedMessage, pureOk] = HardCascade.pure_rs_decode(pureCodec, pureBits);
        assert(pureOk && isequal(decodedMessage, rsMessage), ...
            'HardCascade:PureRsNoiseless', 'Pure-RS baseline failed noiseless decoding.');

        pamBits = [0 0 0 1 1 1 1 0];
        assert(isequal(HardCascade.bits_to_pam4(pamBits), [-3 -1 1 3]) && ...
            isequal(HardCascade.pam4_to_bits_hard([-3 -1 1 3]), pamBits), ...
            'HardCascade:Pam4Mapping', 'PAM4 Gray map/hard demapper mismatch.');
        assert(abs(HardCascade.sigma_from_ebn0_pam4(0, 0.5) - sqrt(5/8)) < 1e-12, ...
            'HardCascade:Pam4Sigma', 'PAM4 sigma convention differs from Python.');

        cfg255 = HardCascade.hard_cascade_config(8, 239);
        cfg127 = HardCascade.hard_cascade_config(7, 113);
        v1_255 = HardCascade.latency_summary_v1(cfg255);
        v1_127 = HardCascade.latency_summary_v1(cfg127);
        v2_255 = HardCascade.latency_summary_v2(cfg255);
        v2_127 = HardCascade.latency_summary_v2(cfg127);
        assert(v1_255.pure_rs == 19 && v1_255.direct_lagrange == 21 && ...
            v1_127.direct_no_share == 20 && v2_255.direct_v2 == 20 && ...
            v2_127.direct_v2 == 17, 'HardCascade:LatencyModel', ...
            'Latency-model reproduction mismatch.');
        [ratio, passes] = HardCascade.kpi_check(v2_255.direct_v2, v2_255.pure_rs);
        assert(passes && ratio <= 10, 'HardCascade:Kpi', 'v2 n255 KPI must pass.');

        report = struct('name', {'bch_conventional', 'bch_direct', ...
            'bch_uncorrectable', 'rs_bm', 'cascade_direct', ...
            'cascade_conventional', 'pure_rs', 'pam4', 'latency'}, ...
            'pass', {decoderPass(1), decoderPass(2), true, true, true, true, true, true, true});
    end
end
end
