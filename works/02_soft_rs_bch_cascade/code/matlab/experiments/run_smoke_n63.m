function all_results = run_smoke_n63(mode)
%RUN_SMOKE_N63 Work 02 n=63 end-to-end reproduction experiment.
%   MODE is 'full' (Python-reference parameters) or 'quick' (CI smoke).
if nargin < 1
    mode = 'full';
end
quick = strcmpi(mode, 'quick');
assert(quick || strcmpi(mode, 'full'), 'run_smoke_n63:Mode', 'mode must be full or quick.');
paths = cascade.ensure_asset_dirs();
cfg = cascade.CascadeConfig(6, 57, 1, 'llosd_tau', 2, 'lcc_eta', 3);
fprintf('[SMOKE TEST] %s\n', cfg.describe());
if quick
    ebn0_list = 8.0;
    minErrors = 99;
    maxFrames = 1;
    suffix = '_matlab_quick';
else
    ebn0_list = 4.0:0.5:8.0;
    minErrors = 30;
    maxFrames = 1500;
    suffix = '_matlab';
end

methods = {'pure_rs_bm', 'Pure RS-BM (hard)'; ...
    'pure_rs_lccbr', 'Pure RS-LCC-BR (soft)'; ...
    'scheme_a', 'Cascade A (LLOSD+LCC-BR)'; ...
    'scheme_b', 'Cascade B (Lagrange shared)'};
all_results = struct();
for ii = 1:size(methods, 1)
    method = methods{ii, 1};
    label = methods{ii, 2};
    fprintf('\n--- %s ---\n', label);
    pack = cascade.get_codec_pack(cfg, method);
    result = cascade.run_bench(label, pack, ebn0_list, cfg.k_rs, cfg.m, ...
        'min_frame_errors', minErrors, 'max_frames', maxFrames, 'verbose', true);
    all_results.(method) = cascade.result_to_struct(result);
end
cascade.write_json(fullfile(paths.data, ['smoke_n63', suffix, '.json']), all_results);

fig = figure('Visible', 'off', 'Position', [100, 100, 700, 500]);
hold on;
styles = {'pure_rs_bm', 'd-', [0.85 0.1 0.1], 'Pure RS-BM'; ...
    'pure_rs_lccbr', 's-', [0.9 0.5 0.05], 'Pure RS-LCC-BR'; ...
    'scheme_a', '^-', [0.1 0.25 0.85], 'Cascade A'; ...
    'scheme_b', 'v--', [0.1 0.55 0.2], 'Cascade B'};
for ii = 1:size(styles, 1)
    result = all_results.(styles{ii, 1});
    semilogy(result.ebn0_db, max(result.fer, 1e-6), styles{ii, 2}, ...
        'Color', styles{ii, 3}, 'DisplayName', styles{ii, 4}, 'MarkerSize', 6);
end
xlabel('Eb/N0 (dB)'); ylabel('FER'); grid on; legend('Location', 'best');
title(sprintf('Smoke test: %s', cfg.describe()));
cascade.save_figure_pair(fig, fullfile(paths.figures, ['smoke_n63', suffix]));
close(fig);
fprintf('Saved data/ and figures/ smoke_n63%s.*\n', suffix);
end
