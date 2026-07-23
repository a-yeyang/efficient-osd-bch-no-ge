function all_results = run_all_configs(mode)
%RUN_ALL_CONFIGS All three paper-reproduction configurations and four methods.
if nargin < 1
    mode = 'full';
end
quick = strcmpi(mode, 'quick');
assert(quick || strcmpi(mode, 'full'), 'run_all_configs:Mode', 'mode must be full or quick.');
paths = cascade.ensure_asset_dirs();
configs = struct();
configs.n255_A = cascade.CascadeConfig(8, 239, 2, 'llosd_tau', 2, 'lcc_eta', 4);
configs.n255_B = cascade.CascadeConfig(8, 235, 1, 'llosd_tau', 1, 'lcc_eta', 4);
configs.n127 = cascade.CascadeConfig(7, 119, 1, 'llosd_tau', 1, 'lcc_eta', 4);
if quick
    minErrors = 99;
    maxFrames = 1;
    suffix = '_matlab_quick';
else
    minErrors = 15;
    maxFrames = 200;
    suffix = '_matlab';
end
methodNames = {'pure_rs_bm', 'pure_rs_lccbr', 'scheme_a', 'scheme_b'};
configNames = fieldnames(configs);
all_results = struct();
for cc = 1:numel(configNames)
    cfgName = configNames{cc};
    cfg = configs.(cfgName);
    fprintf('\n########## %s: %s ##########\n', cfgName, cfg.describe());
    if quick
        if cfg.m == 8
            ebn0_list = 12.0;
        else
            ebn0_list = 11.0;
        end
    elseif cfg.m == 8
        ebn0_list = 6.0:0.5:12.0;
    else
        ebn0_list = 5.0:0.5:11.0;
    end
    perMethod = struct();
    for mm = 1:numel(methodNames)
        method = methodNames{mm};
        pack = cascade.get_codec_pack(cfg, method);
        perMethod.(method) = cascade.run_bench([cfgName, '/', method], pack, ebn0_list, ...
            cfg.k_rs, cfg.m, 'min_frame_errors', minErrors, 'max_frames', maxFrames, 'verbose', true);
    end
    all_results.(cfgName) = perMethod;
    cascade.write_json(fullfile(paths.data, [cfgName, '_all_methods', suffix, '.json']), perMethod);
end

fig = figure('Visible', 'off', 'Position', [100, 100, 1500, 430]);
styles = {'pure_rs_bm', 'd-', [0.85 0.1 0.1], 'Pure RS-BM'; ...
    'pure_rs_lccbr', 's-', [0.9 0.5 0.05], 'Pure RS-LCC-BR'; ...
    'scheme_a', '^-', [0.1 0.25 0.85], 'Cascade A'; ...
    'scheme_b', 'v--', [0.1 0.55 0.2], 'Cascade B'};
for cc = 1:numel(configNames)
    cfgName = configNames{cc};
    cfg = configs.(cfgName);
    subplot(1, numel(configNames), cc); hold on;
    for mm = 1:size(styles, 1)
        result = all_results.(cfgName).(styles{mm, 1});
        semilogy(result.ebn0_db, max(result.fer, 1e-5), styles{mm, 2}, ...
            'Color', styles{mm, 3}, 'DisplayName', styles{mm, 4}, 'MarkerSize', 5);
    end
    xlabel('Eb/N0 (dB)'); ylabel('FER'); grid on;
    title(sprintf('%s: %s', cfgName, cfg.describe()), 'Interpreter', 'none');
    legend('Location', 'best', 'FontSize', 8);
end
cascade.save_figure_pair(fig, fullfile(paths.figures, ['all_configs_fer', suffix]));
close(fig);
cascade.write_json(fullfile(paths.data, ['all_results', suffix, '.json']), all_results);
fprintf('Saved all-configuration MATLAB data and figures.\n');
end
