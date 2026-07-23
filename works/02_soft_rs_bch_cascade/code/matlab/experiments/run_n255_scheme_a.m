function all_results = run_n255_scheme_a(mode)
%RUN_N255_SCHEME_A Main n=255 Scheme-A closed-loop experiment.
%   Full mode matches the Python source; quick mode runs one high-SNR frame.
if nargin < 1
    mode = 'full';
end
quick = strcmpi(mode, 'quick');
assert(quick || strcmpi(mode, 'full'), 'run_n255_scheme_a:Mode', 'mode must be full or quick.');
paths = cascade.ensure_asset_dirs();
cfg = cascade.CascadeConfig(8, 239, 2, 'llosd_tau', 2, 'lcc_eta', 4);
fprintf('[Main experiment: n=255 Scheme A] %s\n', cfg.describe());
if quick
    ebn0_list = 12.0;
    minErrors = 99;
    maxFramesPure = 1;
    maxFramesCascade = 1;
    suffix = '_matlab_quick';
else
    ebn0_list = 6.0:0.5:12.0;
    minErrors = 15;
    maxFramesPure = 400;
    maxFramesCascade = 200;
    suffix = '_matlab';
end
methods = {'pure_rs_bm', 'Pure RS-BM'; ...
    'pure_rs_lccbr', 'Pure RS-LCC-BR'; ...
    'scheme_a', 'Cascade A (LLOSD+LCC-BR)'};
all_results = struct();
for ii = 1:size(methods, 1)
    method = methods{ii, 1};
    label = methods{ii, 2};
    if startsWith(method, 'pure')
        maxFrames = maxFramesPure;
    else
        maxFrames = maxFramesCascade;
    end
    fprintf('\n--- %s ---\n', label);
    pack = cascade.get_codec_pack(cfg, method);
    all_results.(method) = cascade.run_bench(label, pack, ebn0_list, cfg.k_rs, cfg.m, ...
        'min_frame_errors', minErrors, 'max_frames', maxFrames, 'verbose', true);
end
cascade.write_json(fullfile(paths.data, ['n255_scheme_a', suffix, '.json']), all_results);

styles = {'pure_rs_bm', 'd-', [0.85 0.1 0.1], 'Pure RS-BM'; ...
    'pure_rs_lccbr', 's-', [0.9 0.5 0.05], 'Pure RS-LCC-BR'; ...
    'scheme_a', '^-', [0.1 0.25 0.85], 'Cascade A'};
for plotKind = 1:2
    fig = figure('Visible', 'off', 'Position', [100, 100, 700, 500]);
    hold on;
    for ii = 1:size(styles, 1)
        result = all_results.(styles{ii, 1});
        if plotKind == 1
            values = max(result.fer, 1e-5);
            ylabelText = 'FER';
            filename = ['n255_scheme_a_fer', suffix];
        else
            values = max(result.ber, 1e-6);
            ylabelText = 'BER (info bits)';
            filename = ['n255_scheme_a_ber', suffix];
        end
        semilogy(result.ebn0_db, values, styles{ii, 2}, 'Color', styles{ii, 3}, ...
            'DisplayName', styles{ii, 4}, 'MarkerSize', 6);
    end
    xlabel('Eb/N0 (dB)'); ylabel(ylabelText); grid on; legend('Location', 'best');
    title(sprintf('n=255 minimum closed loop: %s', cfg.describe()));
    cascade.save_figure_pair(fig, fullfile(paths.figures, filename));
    close(fig);
end
fprintf('Saved n255 Scheme-A MATLAB data and figures.\n');
end
