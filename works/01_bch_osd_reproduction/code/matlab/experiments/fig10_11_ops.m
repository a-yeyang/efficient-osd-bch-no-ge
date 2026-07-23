function out = fig10_11_ops(options)
%FIG10_11_OPS Reproduce average TEP/TV counts in Fig. 10 and Fig. 11.

if nargin<1, options=struct(); end
out=struct('fig10',fig10(options),'fig11',fig11(options));
end

function out=fig10(options)
paths=work01_setup(); visible=local_opt(options,'visible','off'); scale=local_opt(options,'n_trials_scale',1);
ebn0=3:0.5:7; code31=work01.Core.bch_code(5,2); code63=work01.Core.bch_code(6,3);
r31t1=work01.Experiment.avg_stat(code31,@(c,L)work01.Core.llosd_fast(c,L,1),ebn0,local_n(400,scale),'n_teps',0,'(31,21) tau=1');
r31t2=work01.Experiment.avg_stat(code31,@(c,L)work01.Core.llosd_fast(c,L,2),ebn0,local_n(400,scale),'n_teps',0,'(31,21) tau=2');
r63t1=work01.Experiment.avg_stat(code63,@(c,L)work01.Core.llosd_fast(c,L,1),ebn0,local_n(200,scale),'n_teps',0,'(63,45) tau=1');
r63t2=work01.Experiment.avg_stat(code63,@(c,L)work01.Core.llosd_fast(c,L,2),ebn0,local_n(200,scale),'n_teps',0,'(63,45) tau=2');
fig=figure('Visible',visible,'Position',[100,100,600,400]); hold on;
semilogy(ebn0,r63t2,'s-k','DisplayName','(63,45) BCH, \tau=2'); semilogy(ebn0,r63t1,'o-k','DisplayName','(63,45) BCH, \tau=1');
semilogy(ebn0,r31t2,'s--','Color',[.25 .41 .88],'DisplayName','(31,21) BCH, \tau=2'); semilogy(ebn0,r31t1,'o--','Color',[.25 .41 .88],'DisplayName','(31,21) BCH, \tau=1');
xlabel('Eb/N0 (dB)'); ylabel('Average N_{TEPs}'); grid on; legend('FontSize',8); title('Fig. 10: avg N_{TEPs} processed in LLOSD');
save_figure_pair(fig,fullfile(paths.figures_dir,'fig10_nteps')); close(fig);
out=containers.Map({'ebn0_db','n_teps_63_45_tau1','n_teps_63_45_tau2','n_teps_31_21_tau1','n_teps_31_21_tau2'},{ebn0,r63t1,r63t2,r31t1,r31t2});
write_json_file(fullfile(paths.data_dir,'fig10_nteps.json'),out);
end

function out=fig11(options)
paths=work01_setup(); visible=local_opt(options,'visible','off'); scale=local_opt(options,'n_trials_scale',1);
ebn0=2:0.5:6; code63=work01.Core.bch_code(6,4); code255=work01.Core.bch_code(8,4);
fprintf('BCHCode(6,4) has k = %d\n',code63.k);
r63e4=work01.Experiment.avg_stat(code63,@(c,L)work01.Core.hsd_fast(c,L,1,4),ebn0,local_n(100,scale),'n_tvs',0,'(63,39) eta=4');
r63e6=work01.Experiment.avg_stat(code63,@(c,L)work01.Core.hsd_fast(c,L,1,6),ebn0,local_n(100,scale),'n_tvs',0,'(63,39) eta=6');
r255e4=work01.Experiment.avg_stat(code255,@(c,L)work01.Core.hsd_fast(c,L,1,4),ebn0,local_n(30,scale),'n_tvs',0,'(255,223) eta=4');
r255e6=work01.Experiment.avg_stat(code255,@(c,L)work01.Core.hsd_fast(c,L,1,6),ebn0,local_n(30,scale),'n_tvs',0,'(255,223) eta=6');
fig=figure('Visible',visible,'Position',[100,100,600,400]); hold on;
plot(ebn0,r63e6,'s-k','DisplayName','(63,39) BCH, \eta=6'); plot(ebn0,r63e4,'o-k','DisplayName','(63,39) BCH, \eta=4');
plot(ebn0,r255e6,'s--','Color',[.25 .41 .88],'DisplayName','(255,223) BCH, \eta=6'); plot(ebn0,r255e4,'o--','Color',[.25 .41 .88],'DisplayName','(255,223) BCH, \eta=4');
xlabel('Eb/N0 (dB)'); ylabel('Average N_{TVs}'); grid on; legend('FontSize',8); title('Fig. 11: avg N_{TVs} processed in HSD');
save_figure_pair(fig,fullfile(paths.figures_dir,'fig11_ntvs')); close(fig);
out=containers.Map({'ebn0_db','n_tvs_63_39_eta4','n_tvs_63_39_eta6','n_tvs_255_223_eta4','n_tvs_255_223_eta6'},{ebn0,r63e4,r63e6,r255e4,r255e6});
write_json_file(fullfile(paths.data_dir,'fig11_ntvs.json'),out);
end

function value=local_opt(s,name,default), if isfield(s,name), value=s.(name); else, value=default; end, end
function n=local_n(base,scale), n=max(1,round(base*scale)); end
