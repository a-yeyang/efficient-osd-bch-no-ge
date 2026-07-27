% replot_latency.m — regenerate the latency KPI figure from saved data
% (matlab_results.mat), using the fixed log-scale rendering. The underlying
% numbers are identical to what main() produced; only the figure changes.
here = fileparts(mfilename('fullpath'));
S = load(fullfile(here, 'data', 'matlab_results.mat'));
rbm = S.all_results.n255;

ebn0 = rbm.pure_rs_bm.ebn0_db;
common = min([numel(rbm.pure_rs_bm.avg_f2m_ops), ...
              numel(rbm.pure_rs_lccbr.avg_f2m_ops), ...
              numel(rbm.cascade_llosd.avg_f2m_ops), ...
              numel(rbm.cascade_osd.avg_f2m_ops)]);
casc_fer = rbm.cascade_llosd.fer;
npts = min(common, numel(casc_fer));
cand = find(casc_fer(1:npts) < 0.5);
if isempty(cand), cand = npts; end
key_idx = cand(max(1, numel(cand)-2):end);

lcc  = rbm.pure_rs_lccbr.avg_f2m_ops;
cllo = rbm.cascade_llosd.avg_f2m_ops;
cosd = rbm.cascade_osd.avg_f2m_ops;
cllo_tot = cllo + rbm.cascade_llosd.avg_f2_ops;
cosd_tot = cosd + rbm.cascade_osd.avg_f2_ops;

snrs = ebn0(key_idx);
ra = 100 * (cllo(key_idx) - lcc(key_idx)) ./ max(lcc(key_idx), 1);
ro = 100 * (cosd(key_idx) - lcc(key_idx)) ./ max(lcc(key_idx), 1);

f = figure('Visible', 'off', 'Position', [100 100 900 420]);
subplot(1, 2, 1);
b = bar(categorical(arrayfun(@(s) sprintf('%.1f dB', s), snrs, 'uni', 0)), [ra(:), ro(:)]);
b(1).FaceColor = [0.0 0.6 0.3]; b(2).FaceColor = [0.6 0.6 0.6];
hold on; yline(10, 'k--', '+10% KPI', 'LineWidth', 1.2);
ylabel('时延增幅 vs 纯 RS-LCC-BR (%)');
title('级联时延增幅 (F_2^m 运算量)');
legend({'LLOSD (Lagrange)', 'OSD (Gauss elim)'}, 'Location', 'best', 'FontSize', 8);
grid on;

subplot(1, 2, 2);
b2 = bar(categorical(arrayfun(@(s) sprintf('%.1f dB', s), snrs, 'uni', 0)), ...
    [cllo_tot(key_idx)', cosd_tot(key_idx)']);
b2(1).FaceColor = [0.0 0.3 0.85]; b2(2).FaceColor = [0.6 0.6 0.6];
set(gca, 'YScale', 'log');
ylabel('平均总运算量 / 帧 (F_2^m + F_2, log)');
title({'LLOSD(Lagrange) << OSD(高斯消元)', '(对数轴)'}, 'FontSize', 10);
legend({'LLOSD (Lagrange)', 'OSD (Gauss elim)'}, 'Location', 'northwest', 'FontSize', 8);
grid on;
for j = 1:numel(key_idx)
    idx = key_idx(j);
    spd = cosd_tot(idx) / max(cllo_tot(idx), 1);
    text(j, cosd_tot(idx) * 1.3, sprintf('%.0f\\times', spd), ...
        'HorizontalAlignment', 'center', 'FontSize', 9, ...
        'FontWeight', 'bold', 'Color', [0.0 0.3 0.85]);
end

saveas(f, fullfile(here, 'figures', 'matlab_latency_bars.png'));
savefig(f, fullfile(here, 'figures', 'matlab_latency_bars.fig'));
close(f);
fprintf('REPLOT OK\n');
