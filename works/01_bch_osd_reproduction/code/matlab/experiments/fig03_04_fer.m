function out = fig03_04_fer(options)
%FIG03_04_FER Reproduce the Fig. 3 and Fig. 4 FER curves.
% Options: max_frames_scale (default 1), min_errors_scale (default 1), visible.

if nargin < 1, options = struct(); end
out = struct('fig3',fig3_31_21(options),'fig4',fig4_63_45(options));
end

function out = fig3_31_21(options)
paths = work01_setup(); scale = local_opt(options,'max_frames_scale',1);
minScale = local_opt(options,'min_errors_scale',1); visible = local_opt(options,'visible','off');
code = work01.Core.bch_code(5,2); ebn0 = 2:1:8; me = max(1,round(60*minScale));
fprintf('\n=== Fig 3: (31, 21) BCH ===\n');
results = containers.Map;
results('BM') = work01.Experiment.run_fer(code,@work01.Core.bm_wrap,ebn0,me,max(1,round(20000*scale)),0,'BM');
results('OSD(1)') = work01.Experiment.run_fer(code,@(c,L) work01.Core.osd_decode(c,L,1),ebn0,me,max(1,round(15000*scale)),0,'OSD(1)');
results('LLOSD(1)') = work01.Experiment.run_fer(code,@(c,L) work01.Core.llosd_fast(c,L,1),ebn0,me,max(1,round(15000*scale)),0,'LLOSD(1)');
results('LLOSD(2)') = work01.Experiment.run_fer(code,@(c,L) work01.Core.llosd_fast(c,L,2),ebn0,me,max(1,round(15000*scale)),0,'LLOSD(2)');
results('ML') = work01.Experiment.run_fer(code,@(c,L) work01.Core.llosd_fast(c,L,4),ebn0,me,max(1,round(15000*scale)),0,'ML(~LLOSD 4)');
fig = figure('Visible',visible,'Position',[100,100,600,450]); hold on;
styles = {'BM','d-','r';'OSD(1)','o-','k';'LLOSD(1)','v-','k';'LLOSD(2)','^-','k';'ML','--','r'};
for i=1:size(styles,1), semilogy(ebn0,max(results(styles{i,1}),1e-6),styles{i,2},'Color',styles{i,3},'DisplayName',styles{i,1},'MarkerSize',5); end
xlabel('Eb/N0 (dB)'); ylabel('FER'); grid on; legend('FontSize',9); title('Fig. 3: (31, 21) BCH code');
save_figure_pair(fig,fullfile(paths.figures_dir,'fig03_31_21')); close(fig);
out = containers.Map([{'ebn0_db'},results.keys],[{ebn0},results.values]);
write_json_file(fullfile(paths.data_dir,'fig03_31_21.json'),out);
end

function out = fig4_63_45(options)
paths = work01_setup(); scale = local_opt(options,'max_frames_scale',1);
minScale = local_opt(options,'min_errors_scale',1); visible = local_opt(options,'visible','off');
code = work01.Core.bch_code(6,3); ebn0 = 1:1:7; me = max(1,round(60*minScale));
fprintf('\n=== Fig 4: (63, 45) BCH ===\n'); results = containers.Map;
results('BM') = work01.Experiment.run_fer(code,@work01.Core.bm_wrap,ebn0,me,max(1,round(20000*scale)),0,'BM');
results('OSD(1)') = work01.Experiment.run_fer(code,@(c,L) work01.Core.osd_decode(c,L,1),ebn0,me,max(1,round(15000*scale)),0,'OSD(1)');
results('LLOSD(1)') = work01.Experiment.run_fer(code,@(c,L) work01.Core.llosd_fast(c,L,1),ebn0,me,max(1,round(15000*scale)),0,'LLOSD(1)');
results('LLOSD(2)') = work01.Experiment.run_fer(code,@(c,L) work01.Core.llosd_fast(c,L,2),ebn0,me,max(1,round(15000*scale)),0,'LLOSD(2)');
results('LLOSD(3)') = work01.Experiment.run_fer(code,@(c,L) work01.Core.llosd_fast(c,L,3),ebn0,me,max(1,round(15000*scale)),0,'LLOSD(3)');
results('ML') = work01.Experiment.run_fer(code,@(c,L) work01.Core.llosd_fast(c,L,4),ebn0,me,max(1,round(10000*scale)),0,'ML(~LLOSD 4)');
fig = figure('Visible',visible,'Position',[100,100,600,450]); hold on;
styles = {'BM','d-','r';'LLOSD(1)','v-','k';'LLOSD(2)','^-','k';'LLOSD(3)','x-','k';'OSD(1)','o-','r';'ML','--','r'};
for i=1:size(styles,1), semilogy(ebn0,max(results(styles{i,1}),1e-7),styles{i,2},'Color',styles{i,3},'DisplayName',styles{i,1},'MarkerSize',5); end
xlabel('Eb/N0 (dB)'); ylabel('FER'); grid on; legend('FontSize',9); title('Fig. 4: (63, 45) BCH code');
save_figure_pair(fig,fullfile(paths.figures_dir,'fig04_63_45')); close(fig);
out = containers.Map([{'ebn0_db'},results.keys],[{ebn0},results.values]);
write_json_file(fullfile(paths.data_dir,'fig04_63_45.json'),out);
end

function value = local_opt(s,name,default)
if isfield(s,name), value=s.(name); else, value=default; end
end
