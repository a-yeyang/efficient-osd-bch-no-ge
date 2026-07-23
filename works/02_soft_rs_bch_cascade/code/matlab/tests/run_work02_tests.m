function run_work02_tests()
%RUN_WORK02_TESTS Fast, deterministic assertion suite for the MATLAB port.
%   This test intentionally uses short, high-confidence paths where exact
%   decoded outputs are known, plus a tiny noisy Monte Carlo smoke test.
thisFile = mfilename('fullpath');
matlabRoot = fileparts(fileparts(thisFile));
addpath(matlabRoot);
addpath(fullfile(matlabRoot, 'experiments'));
set(groot, 'defaultFigureVisible', 'off');
rng(20260723, 'twister');

requiredEntrypoints = {'run_smoke_n63', 'run_n255_scheme_a', 'run_all_configs', ...
    'run_scheme_ab_savings', 'run_kpi_analysis'};
for ii = 1:numel(requiredEntrypoints)
    assert(exist(requiredEntrypoints{ii}, 'file') == 2, ...
        'Missing required experiment entry point: %s', requiredEntrypoints{ii});
end

% GF(2^m) arithmetic and BCH construction / bounded-distance decoding.
gf = cascade.GF(4);
assert(gf.mul(3, 5) == 15, 'GF multiplication table mismatch.');
assert(gf.div(gf.mul(7, 13), 13) == 7, 'GF division inverse mismatch.');
assert(gf.poly_eval([1, 1, 1], 2) == cascade.gfxor(cascade.gfxor(1, 2), gf.mul(2, 2)), ...
    'GF polynomial evaluation mismatch.');
bch = cascade.BCHCode(6, 1);
bchMessage = randi([0, 1], 1, bch.k);
bchCodeword = bch.encode(bchMessage);
receivedBch = bchCodeword;
receivedBch(11) = 1 - receivedBch(11);
[decodedBch, bchOK] = bch.bm_decode(receivedBch);
assert(bchOK && isequal(decodedBch, bchCodeword), 'BCH BM failed a correctable one-bit error.');
[llosdWord, llosdStats] = cascade.llosd_fast(bch, (1 - 2 * bchCodeword) * 30, 2, true, true);
assert(isequal(llosdWord, bchCodeword), 'LLOSD did not recover a clean BCH frame.');
assert(llosdStats.n_teps >= 1 && llosdStats.n_bch_candidates >= 1, 'LLOSD statistics were not populated.');

% Gray PAM4 map/hard decision/LLR sign convention and AWGN sigma formula.
pamBits = [0 0 0 1 1 1 1 0];
pamSymbols = cascade.bits_to_pam4(pamBits);
assert(isequal(pamSymbols, [-3 -1 1 3]), 'PAM4 Gray mapping mismatch.');
assert(isequal(cascade.pam4_to_bits_hard(pamSymbols), pamBits), 'PAM4 hard demapping mismatch.');
llr = cascade.pam4_bit_llr(pamSymbols, 0.2);
assert(all(sign(llr) == [1 1 1 -1 -1 -1 -1 1]), 'PAM4 LLR convention mismatch.');
assert(abs(cascade.sigma_from_ebn0_pam4(0, 1) - sqrt(5 / 4)) < 1e-12, 'PAM4 sigma mismatch.');

% RS systematic encoding, BM correction, and LCC-BR Chase path.
rs = cascade.RSCode(6, 57);
rsMessage = randi([0, 63], 1, rs.k);
rsCodeword = rs.encode_systematic(rsMessage);
receivedRs = rsCodeword;
errorPositions = [2, 18, 42];
receivedRs(errorPositions) = cascade.gfxor(receivedRs(errorPositions), [1, 7, 21]);
[decodedRs, rsOK] = rs.bm_decode(receivedRs);
assert(rsOK && isequal(decodedRs, rsCodeword), 'RS BM failed t correctable symbol errors.');
reliability = ones(1, rs.n) * 10;
reliability(errorPositions) = 0.01;
[softRs, softOK] = rs.lcc_br_decode(receivedRs, reliability, 3);
assert(softOK && isequal(softRs, rsCodeword), 'RS LCC-BR failed the Chase correction case.');

% Lagrange cache basis identities and real cache reuse accounting.
cache = cascade.LagrangeCache(gf, gf.n);
support = [1, 2, 4];
assert(cache.lagrange_basis(1, support, 1) == 1, 'Lagrange basis is not one at its support point.');
assert(cache.lagrange_basis(1, support, 2) == 0, 'Lagrange basis is not zero at another support point.');
cache.denominator_product(1, support);
cache.denominator_product(1, support);
cacheStats = cache.stats();
assert(cacheStats.cache_size >= 1 && cacheStats.ops_saved >= 2, 'Lagrange cache did not record reuse.');
assert(isequal(size(cache.build_pairwise_diff()), [gf.n, gf.n]), 'Lagrange pairwise table shape mismatch.');

% Complete Scheme A / B recovery and reported cache-operation savings.
cfg = cascade.CascadeConfig(6, 57, 1, 'llosd_tau', 2, 'lcc_eta', 3);
cascadeCodec = cascade.CascadedCodec(cfg);
cascadeMessage = randi([0, 63], 1, cfg.k_rs);
cascadeBits = cascadeCodec.encode(cascadeMessage);
cleanLlr = (1 - 2 * cascadeBits) * 30;
[decodedA, resultA] = cascadeCodec.decode(cleanLlr, 'scheme_a');
[decodedB, resultB] = cascadeCodec.decode(cleanLlr, 'scheme_b');
assert(isequal(decodedA, cascadeMessage), 'Scheme A did not recover a clean cascade frame.');
assert(isequal(decodedB, cascadeMessage), 'Scheme B did not recover a clean cascade frame.');
innerGeneratorCost = 2 * (cascadeCodec.bch.n^2 - ...
    (cascadeCodec.bch.n - 2 * cascadeCodec.bch.t)^2 + ...
    (cascadeCodec.bch.n - 2 * cascadeCodec.bch.t));
assert(resultA.counters.f2m >= cascadeCodec.n_bch_blocks * innerGeneratorCost, ...
    'Scheme A operation count omitted one or more inner LLOSD generator costs.');
assert(resultB.counters.f2m <= resultA.counters.f2m, 'Scheme B did not report cache savings.');
assert(isfield(resultB, 'cache_stats'), 'Scheme B result lacks cache statistics.');
assert(resultB.cache_stats.cache_size > 0 && resultB.cache_stats.ops_saved > 0, ...
    'Scheme B did not create and reuse Lagrange denominator-cache entries.');
assert(resultB.cache_stats.basis_queries > 0 && resultB.cache_stats.alpha_hits > 0 && ...
    resultB.cache_stats.pairwise_table_built, ...
    'Scheme B did not use the shared cache in both Lagrange and RS algebra paths.');
assert(resultB.cache_accounted_savings > 0 && ...
    resultA.counters.f2m - resultB.counters.f2m == resultB.cache_accounted_savings, ...
    'Scheme B accounting is not tied to observed cache reuse.');

% Tiny noisy Monte Carlo run ensures the simulator and function-handle pack work.
pack = cascade.get_codec_pack(cfg, 'pure_rs_bm');
benchmark = cascade.run_bench('test/pure_rs_bm', pack, 8, cfg.k_rs, cfg.m, ...
    'seed', 7, 'min_frame_errors', 99, 'max_frames', 2, 'verbose', false);
assert(benchmark.n_frames == 2 && numel(benchmark.fer) == 1, 'Monte Carlo runner did not complete requested frames.');
assert(isfinite(benchmark.avg_f2m_ops) && isfinite(benchmark.avg_latency_us), 'Monte Carlo metrics are invalid.');

% A quick public experiment validates path independence and asset generation.
smokeResults = run_smoke_n63('quick');
assert(isfield(smokeResults, 'scheme_a') && isfield(smokeResults, 'scheme_b'), 'Smoke experiment omitted cascade methods.');
paths = cascade.work02_paths();
assert(isfile(fullfile(paths.data, 'smoke_n63_matlab_quick.json')), 'Smoke experiment did not write its JSON result.');
assert(isfile(fullfile(paths.figures, 'smoke_n63_matlab_quick.png')), 'Smoke experiment did not write its PNG figure.');

cascade.write_json(fullfile(paths.data, 'matlab_work02_test_result.json'), ...
    struct('passed', true, 'timestamp', datestr(now, 30)));
fprintf('Work 02 MATLAB tests passed.\n');
end
