function savings_by_config = run_scheme_ab_savings(mode)
%RUN_SCHEME_AB_SAVINGS Plot GF-operation savings from Scheme-B cache sharing.
if nargin < 1
    mode = 'full';
end
quick = strcmpi(mode, 'quick');
assert(quick || strcmpi(mode, 'full'), 'run_scheme_ab_savings:Mode', 'mode must be full or quick.');
paths = cascade.ensure_asset_dirs();
if quick
    suffix = '_matlab_quick';
else
    suffix = '_matlab';
end
source = fullfile(paths.data, ['all_results', suffix, '.json']);
if ~isfile(source)
    run_all_configs(mode);
end
data = cascade.read_json(source);
configNames = {'n255_A', 'n255_B', 'n127'};
savings_by_config = struct();
fig = figure('Visible', 'off', 'Position', [100, 100, 1500, 440]);
for ii = 1:numel(configNames)
    name = configNames{ii};
    result = data.(name);
    opsA = result.scheme_a.avg_f2m_ops;
    opsB = result.scheme_b.avg_f2m_ops;
    savings = zeros(size(opsA));
    nonzero = opsA > 0;
    savings(nonzero) = (opsA(nonzero) - opsB(nonzero)) ./ opsA(nonzero) * 100;
    savings_by_config.(name) = struct('ebn0_db', result.scheme_a.ebn0_db, ...
        'ops_a', opsA, 'ops_b', opsB, 'savings_percent', savings);
    fprintf('\n%s:\n', name);
    for jj = 1:numel(savings)
        fprintf('  @ %.1f dB: A=%.0f, B=%.0f, saved=%.2f%%\n', ...
            result.scheme_a.ebn0_db(jj), opsA(jj), opsB(jj), savings(jj));
    end
    subplot(1, numel(configNames), ii);
    bar(1:numel(savings), savings, 'FaceColor', [0.1 0.55 0.2], 'FaceAlpha', 0.7);
    xticks(1:numel(savings));
    xticklabels(compose('%.1f', result.scheme_a.ebn0_db));
    xlabel('Eb/N0 (dB)'); ylabel('Ops savings (%)'); grid on;
    title(sprintf('%s: Scheme B vs A', name), 'Interpreter', 'none');
    yline(0, 'k-');
end
cascade.save_figure_pair(fig, fullfile(paths.figures, ['scheme_ab_savings', suffix]));
close(fig);
cascade.write_json(fullfile(paths.data, ['scheme_ab_savings', suffix, '.json']), savings_by_config);
fprintf('Saved Scheme A/B savings data and figure.\n');
end
