function out = fig05_comparison(options)
%FIG05_COMPARISON Reproduce Fig. 5 with the original simplified baselines.

if nargin < 1, options = struct(); end
paths = work01_setup(); scale=local_opt(options,'max_frames_scale',1);
me=max(1,round(60*local_opt(options,'min_errors_scale',1))); visible=local_opt(options,'visible','off');
code=work01.Core.bch_code(6,3); ebn0=1:1:7; results=containers.Map;
fprintf('=== Fig 5: (63,45) comparison ===\n');
results('LLOSD(3)')=work01.Experiment.run_fer(code,@(c,L) work01.Core.llosd_fast(c,L,3),ebn0,me,max(1,round(15000*scale)),0,'LLOSD(3)');
results('SLLOSD(3,2)')=work01.Experiment.run_fer(code,@(c,L) work01.Core.sllosd_fast(c,L,[3,2]),ebn0,me,max(1,round(15000*scale)),0,'SLLOSD(3,2)');
results('YSVL OSD(1)')=work01.Experiment.run_fer(code,@(c,L) work01.Core.ysvl_osd_decode(c,L,1),ebn0,me,max(1,round(15000*scale)),0,'YSVL(1)');
results('CJ OSD(1)')=work01.Experiment.run_fer(code,@(c,L) work01.Core.cj_osd_decode(c,L,1),ebn0,me,max(1,round(15000*scale)),0,'CJ(1)');
results('OSD(1)')=work01.Experiment.run_fer(code,@(c,L) work01.Core.osd_decode(c,L,1),ebn0,me,max(1,round(15000*scale)),0,'OSD(1)');
results('ML')=work01.Experiment.run_fer(code,@(c,L) work01.Core.llosd_fast(c,L,4),ebn0,me,max(1,round(10000*scale)),0,'ML(~LLOSD 4)');
fig=figure('Visible',visible,'Position',[100,100,700,500]); hold on;
styles={'LLOSD(3)','k^-';'SLLOSD(3,2)','s-';'YSVL OSD(1)','d-';'CJ OSD(1)','x-';'OSD(1)','o-';'ML','r--'};
colors={'k',[0.25,0.41,0.88],[0.85,0.65,0],[1,0.5,0],'r','r'};
for i=1:size(styles,1), semilogy(ebn0,max(results(styles{i,1}),1e-7),styles{i,2},'Color',colors{i},'DisplayName',styles{i,1},'MarkerSize',5); end
xlabel('Eb/N0 (dB)'); ylabel('FER'); grid on; legend('FontSize',9,'Location','southwest'); title('Fig. 5: (63,45) BCH comparison');
save_figure_pair(fig,fullfile(paths.figures_dir,'fig05_63_45_comparison')); close(fig);
out=containers.Map([{'ebn0_db'},results.keys],[{ebn0},results.values]);
write_json_file(fullfile(paths.data_dir,'fig05_63_45.json'),out);
end

function value=local_opt(s,name,default)
if isfield(s,name), value=s.(name); else, value=default; end
end
