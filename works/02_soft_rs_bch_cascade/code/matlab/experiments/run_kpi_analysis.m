function summary = run_kpi_analysis(mode)
%RUN_KPI_ANALYSIS KPI/parallelism analysis for the n=255 Scheme-A study.
if nargin < 1
    mode = 'full';
end
quick = strcmpi(mode, 'quick');
assert(quick || strcmpi(mode, 'full'), 'run_kpi_analysis:Mode', 'mode must be full or quick.');
paths = cascade.ensure_asset_dirs();
if quick
    suffix = '_matlab_quick';
else
    suffix = '_matlab';
end
source = fullfile(paths.data, ['n255_scheme_a', suffix, '.json']);
if ~isfile(source)
    run_n255_scheme_a(mode);
end
data = cascade.read_json(source);
ebn0 = data.pure_rs_bm.ebn0_db;
nPoints = min([numel(ebn0), numel(data.pure_rs_lccbr.ebn0_db), numel(data.scheme_a.ebn0_db)]);
summary = struct('ebn0_db', [], 'kpi_case_a', [], 'kpi_case_b', [], ...
    'rs_bm_ops', [], 'rs_lccbr_ops', [], 'cascade_a_ops', []);
for ii = 1:nPoints
    bm = data.pure_rs_bm.avg_f2m_ops(ii);
    lcc = data.pure_rs_lccbr.avg_f2m_ops(ii);
    cas = data.scheme_a.avg_f2m_ops(ii);
    if bm > 0
        ratioA = (cas - bm) / bm;
    else
        ratioA = inf;
    end
    if lcc > 0
        ratioB = (cas - lcc) / lcc;
    else
        ratioB = inf;
    end
    summary.ebn0_db(end + 1) = ebn0(ii);
    summary.kpi_case_a(end + 1) = ratioA;
    summary.kpi_case_b(end + 1) = ratioB;
    summary.rs_bm_ops(end + 1) = bm;
    summary.rs_lccbr_ops(end + 1) = lcc;
    summary.cascade_a_ops(end + 1) = cas;
    fprintf('@ %.1f dB: BM=%.0f, LCC-BR=%.0f, Cascade A=%.0f; vs BM=%+.1f%%, vs LCC=%+.1f%%\n', ...
        ebn0(ii), bm, lcc, cas, 100 * ratioA, 100 * ratioB);
    for parallelism = [1, 4, 16, 64]
        fprintf('  P=%d: RS-BM=%.0f, RS-LCC-BR=%.0f, Cascade=%.0f cycles\n', ...
            parallelism, bm / parallelism, lcc / parallelism, cas / parallelism);
    end
end
cascade.write_json(fullfile(paths.data, ['kpi_summary', suffix, '.json']), summary);

fig = figure('Visible', 'off', 'Position', [100, 100, 1100, 430]);
subplot(1, 2, 1);
bar(summary.kpi_case_a * 100, 'FaceColor', [0.1 0.55 0.2]); hold on; yline(10, 'k--', '+10%');
xticks(1:numel(summary.ebn0_db)); xticklabels(compose('%.1f', summary.ebn0_db));
xlabel('Eb/N0 (dB)'); ylabel('Latency increase (%)'); title('Case (a): Cascade vs RS-BM'); grid on;
subplot(1, 2, 2);
bar(summary.kpi_case_b * 100, 'FaceColor', [0.1 0.25 0.85]); hold on; yline(10, 'k--', '+10%');
xticks(1:numel(summary.ebn0_db)); xticklabels(compose('%.1f', summary.ebn0_db));
xlabel('Eb/N0 (dB)'); ylabel('Latency increase (%)'); title('Case (b): Cascade vs RS-LCC-BR'); grid on;
cascade.save_figure_pair(fig, fullfile(paths.figures, ['n255_kpi_bars', suffix]));
close(fig);
fprintf('Saved KPI summary and figure.\n');
end
