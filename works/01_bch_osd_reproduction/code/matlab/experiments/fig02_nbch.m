function out = fig02_nbch(options)
%FIG02_NBCH Reproduce Fig. 2: average BCH candidates versus Eb/N0.
% Call with no argument for the Python script's full settings.  An options
% struct may override ebn0_db, n_trials_63, n_trials_31, and visible.

if nargin < 1, options = struct(); end
paths = work01_setup();
ebn0 = local_opt(options,'ebn0_db',2:1:10);
code63 = work01.Core.bch_code(6,3);
code31 = work01.Core.bch_code(5,2);
nbch63 = work01.Experiment.avg_nbch(code63,3,ebn0,local_opt(options,'n_trials_63',200),0);
nbch31 = work01.Experiment.avg_nbch(code31,2,ebn0,local_opt(options,'n_trials_31',800),0);

fig = figure('Visible',local_opt(options,'visible','off'),'Position',[100,100,600,400]);
plot(ebn0,nbch63,'k^-','DisplayName','Simulation, LLOSD (3)'); hold on;
yline(5,'r-','DisplayName','Theoretical, LLOSD (3)');
plot(ebn0,nbch31,'kv--','DisplayName','Simulation, LLOSD (2)');
yline(3,'r--','DisplayName','Theoretical, LLOSD (2)');
text(6,5.7,'(63, 45) BCH','FontSize',9); text(6,3.3,'(31, 21) BCH','FontSize',9);
xlabel('Eb/N0 (dB)'); ylabel('N_{BCH}'); ylim([0,10]); grid on;
legend('FontSize',8,'Location','northeast');
save_figure_pair(fig,fullfile(paths.figures_dir,'fig02_nbch'));
close(fig);

out = containers.Map({'ebn0_db','nbch_63_45_tau3','nbch_31_21_tau2', ...
    'theoretical_63_45','theoretical_31_21'}, {ebn0,nbch63,nbch31,5,3});
write_json_file(fullfile(paths.data_dir,'fig02_nbch.json'),out);
end

function value = local_opt(s,name,default)
if isfield(s,name), value = s.(name); else, value = default; end
end
