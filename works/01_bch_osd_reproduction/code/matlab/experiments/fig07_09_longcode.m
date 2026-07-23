function out = fig07_09_longcode(options)
%FIG07_09_LONGCODE Reproduce Fig. 7 and Fig. 9 long-code FER experiments.
% The no-argument path keeps the Python frame/error budgets.  Use
% max_frames_scale < 1 for a short exploratory run.

if nargin < 1, options = struct(); end
out = struct('fig7',fig7_127_99(options),'fig9',fig9_255_223(options));
end

function out = fig7_127_99(options)
paths=work01_setup(); scale=local_opt(options,'max_frames_scale',1); visible=local_opt(options,'visible','off');
minErrors=max(1,round(40*local_opt(options,'min_errors_scale',1)));
code=work01.Core.bch_code(7,4); ebn0=3:1:7; r=containers.Map;
fprintf('\n=== Fig 7: (%d,%d) ===\n',code.n,code.k);
r('OSD(1)')=work01.Experiment.run_fer(code,@(c,L) work01.Core.osd_decode(c,L,1),ebn0,minErrors,local_frames(3000,scale),0,'OSD(1)');
r('OSD(2)')=work01.Experiment.run_fer(code,@(c,L) work01.Core.osd_decode(c,L,2),ebn0,minErrors,local_frames(3000,scale),0,'OSD(2)');
r('LLOSD(2)')=work01.Experiment.run_fer(code,@(c,L) work01.Core.llosd_fast(c,L,2),ebn0,minErrors,local_frames(3000,scale),0,'LLOSD(2)');
r('LLOSD(3)')=work01.Experiment.run_fer(code,@(c,L) work01.Core.llosd_fast(c,L,3),ebn0,minErrors,local_frames(2000,scale),0,'LLOSD(3)');
r('HSD(1,4)')=work01.Experiment.run_fer(code,@(c,L) work01.Core.hsd_fast(c,L,1,4),ebn0,minErrors,local_frames(3000,scale),0,'HSD(1,4)');
r('HSD(1,6)')=work01.Experiment.run_fer(code,@(c,L) work01.Core.hsd_fast(c,L,1,6),ebn0,minErrors,local_frames(2000,scale),0,'HSD(1,6)');
r('HSD(1,8)')=work01.Experiment.run_fer(code,@(c,L) work01.Core.hsd_fast(c,L,1,8),ebn0,minErrors,local_frames(1500,scale),0,'HSD(1,8)');
fig=figure('Visible',visible,'Position',[100,100,600,450]); hold on;
styles={'OSD(1)','o-','r';'OSD(2)','s-','r';'LLOSD(2)','d-',[.25 .41 .88];'LLOSD(3)','d--',[.25 .41 .88];'HSD(1,4)','^-','k';'HSD(1,6)','v-','k';'HSD(1,8)','x-','k'};
for i=1:size(styles,1), semilogy(ebn0,max(r(styles{i,1}),1e-6),styles{i,2},'Color',styles{i,3},'DisplayName',styles{i,1},'MarkerSize',4); end
xlabel('Eb/N0 (dB)'); ylabel('FER'); grid on; legend('FontSize',8); title('Fig. 7: (127,99) BCH');
save_figure_pair(fig,fullfile(paths.figures_dir,'fig07_127_99')); close(fig);
out=containers.Map([{'ebn0_db'},r.keys],[{ebn0},r.values]); write_json_file(fullfile(paths.data_dir,'fig07_127_99.json'),out);
end

function out = fig9_255_223(options)
paths=work01_setup(); scale=local_opt(options,'max_frames_scale',1); visible=local_opt(options,'visible','off');
minErrors=max(1,round(30*local_opt(options,'min_errors_scale',1)));
code=work01.Core.bch_code(8,4); ebn0=3.5:0.5:7; r=containers.Map;
fprintf('\n=== Fig 9: (%d,%d) ===\n',code.n,code.k);
r('BM')=work01.Experiment.run_fer(code,@work01.Core.bm_wrap,ebn0,minErrors,local_frames(3000,scale),0,'BM');
r('OSD(1)')=work01.Experiment.run_fer(code,@(c,L) work01.Core.osd_decode(c,L,1),ebn0,minErrors,local_frames(1000,scale),0,'OSD(1)');
r('LLOSD(2)')=work01.Experiment.run_fer(code,@(c,L) work01.Core.llosd_fast(c,L,2),ebn0,minErrors,local_frames(800,scale),0,'LLOSD(2)');
r('PLCC(6)')=work01.Experiment.run_fer(code,@(c,L) work01.Core.plcc_decode(c,L,6),ebn0,minErrors,local_frames(500,scale),0,'PLCC(6)');
r('PLCC(8)')=work01.Experiment.run_fer(code,@(c,L) work01.Core.plcc_decode(c,L,8),ebn0,minErrors,local_frames(400,scale),0,'PLCC(8)');
r('HSD(1,4)')=work01.Experiment.run_fer(code,@(c,L) work01.Core.hsd_fast(c,L,1,4),ebn0,minErrors,local_frames(500,scale),0,'HSD(1,4)');
r('HSD(1,6)')=work01.Experiment.run_fer(code,@(c,L) work01.Core.hsd_fast(c,L,1,6),ebn0,minErrors,local_frames(400,scale),0,'HSD(1,6)');
r('HSD(1,8)')=work01.Experiment.run_fer(code,@(c,L) work01.Core.hsd_fast(c,L,1,8),ebn0,minErrors,local_frames(300,scale),0,'HSD(1,8)');
r('ML')=r('HSD(1,8)');
fig=figure('Visible',visible,'Position',[100,100,600,450]); hold on;
styles={'BM','d-','r';'OSD(1)','o-','r';'LLOSD(2)','x-',[.25 .41 .88];'PLCC(6)','s-','k';'PLCC(8)','s--','k';'HSD(1,4)','^-','k';'HSD(1,6)','v-','k';'HSD(1,8)','<-','k';'ML','r--','r'};
for i=1:size(styles,1), semilogy(ebn0,max(r(styles{i,1}),1e-5),styles{i,2},'Color',styles{i,3},'DisplayName',styles{i,1},'MarkerSize',4); end
xlabel('Eb/N0 (dB)'); ylabel('FER'); grid on; legend('FontSize',8,'Location','southwest'); title('Fig. 9: (255,223) BCH');
save_figure_pair(fig,fullfile(paths.figures_dir,'fig09_255_223')); close(fig);
out=containers.Map([{'ebn0_db'},r.keys],[{ebn0},r.values]); write_json_file(fullfile(paths.data_dir,'fig09_255_223.json'),out);
end

function value=local_opt(s,name,default)
if isfield(s,name), value=s.(name); else, value=default; end
end
function n=local_frames(base,scale), n=max(1,round(base*scale)); end
