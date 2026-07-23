function out = fig07_09_hsd_rerun(options)
%FIG07_09_HSD_RERUN Recompute only HSD curves, preserving prior JSON curves.
% Run fig07_09_longcode first; this function intentionally errors if its data
% files are missing, matching the dependency in the Python script.

if nargin < 1, options=struct(); end
out=struct('fig7',rerun7(options),'fig9',rerun9(options));
end

function output=rerun7(options)
paths=work01_setup(); old=jsondecode(fileread(fullfile(paths.data_dir,'fig07_127_99.json')));
scale=local_opt(options,'max_frames_scale',1); visible=local_opt(options,'visible','off'); me=max(1,round(40*local_opt(options,'min_errors_scale',1)));
code=work01.Core.bch_code(7,4); ebn0=3:1:7;
output=local_to_map(old); output('ebn0_db')=ebn0;
output('HSD(1,4)')=work01.Experiment.run_fer(code,@(c,L)work01.Core.hsd_fast(c,L,1,4),ebn0,me,local_frames(3000,scale),0,'HSD(1,4)');
output('HSD(1,6)')=work01.Experiment.run_fer(code,@(c,L)work01.Core.hsd_fast(c,L,1,6),ebn0,me,local_frames(2000,scale),0,'HSD(1,6)');
output('HSD(1,8)')=work01.Experiment.run_fer(code,@(c,L)work01.Core.hsd_fast(c,L,1,8),ebn0,me,local_frames(1500,scale),0,'HSD(1,8)');
local_plot7(output,ebn0,visible,fullfile(paths.figures_dir,'fig07_127_99')); write_json_file(fullfile(paths.data_dir,'fig07_127_99.json'),output);
end

function output=rerun9(options)
paths=work01_setup(); old=jsondecode(fileread(fullfile(paths.data_dir,'fig09_255_223.json')));
scale=local_opt(options,'max_frames_scale',1); visible=local_opt(options,'visible','off'); me=max(1,round(30*local_opt(options,'min_errors_scale',1)));
code=work01.Core.bch_code(8,4); ebn0=3.5:0.5:7;
output=local_to_map(old); output('ebn0_db')=ebn0;
output('HSD(1,4)')=work01.Experiment.run_fer(code,@(c,L)work01.Core.hsd_fast(c,L,1,4),ebn0,me,local_frames(500,scale),0,'HSD(1,4)');
output('HSD(1,6)')=work01.Experiment.run_fer(code,@(c,L)work01.Core.hsd_fast(c,L,1,6),ebn0,me,local_frames(400,scale),0,'HSD(1,6)');
output('HSD(1,8)')=work01.Experiment.run_fer(code,@(c,L)work01.Core.hsd_fast(c,L,1,8),ebn0,me,local_frames(300,scale),0,'HSD(1,8)');
local_plot9(output,ebn0,visible,fullfile(paths.figures_dir,'fig09_255_223')); write_json_file(fullfile(paths.data_dir,'fig09_255_223.json'),output);
end

function m=local_to_map(s)
f=fieldnames(s); m=containers.Map;
% jsondecode turns invalid JSON object keys into makeValidName keys. Restore
% the known published labels before writing the rerun JSON back out.
labels={'ebn0_db','BM','OSD(1)','OSD(2)','LLOSD(2)','LLOSD(3)', ...
    'PLCC(6)','PLCC(8)','HSD(1,4)','HSD(1,6)','HSD(1,8)','ML'};
for i=1:numel(f)
    key=f{i}; restored=key;
    for j=1:numel(labels)
        if strcmp(key,matlab.lang.makeValidName(labels{j}))
            restored=labels{j}; break;
        end
    end
    m(restored)=s.(key);
end
end
function local_plot7(r,x,visible,outbase)
fig=figure('Visible',visible); hold on; names={'OSD_1_','OSD_2_','LLOSD_2_','LLOSD_3_','HSD_1_4_','HSD_1_6_','HSD_1_8_'}; labels={'OSD(1)','OSD(2)','LLOSD(2)','LLOSD(3)','HSD(1,4)','HSD(1,6)','HSD(1,8)'};
for i=1:numel(names), if isKey(r,labels{i}), y=r(labels{i}); elseif isKey(r,names{i}), y=r(names{i}); else, continue; end, semilogy(x,max(y,1e-6),'-o','DisplayName',labels{i}); end
xlabel('Eb/N0 (dB)'); ylabel('FER'); grid on; legend('FontSize',8); title('Fig. 7: fixed HSD'); save_figure_pair(fig,outbase); close(fig);
end
function local_plot9(r,x,visible,outbase)
fig=figure('Visible',visible); hold on; names={'BM','OSD_1_','LLOSD_2_','PLCC_6_','PLCC_8_','HSD_1_4_','HSD_1_6_','HSD_1_8_'}; labels={'BM','OSD(1)','LLOSD(2)','PLCC(6)','PLCC(8)','HSD(1,4)','HSD(1,6)','HSD(1,8)'};
for i=1:numel(names), if isKey(r,labels{i}), y=r(labels{i}); elseif isKey(r,names{i}), y=r(names{i}); else, continue; end, semilogy(x,max(y,1e-5),'-o','DisplayName',labels{i}); end
xlabel('Eb/N0 (dB)'); ylabel('FER'); grid on; legend('FontSize',8); title('Fig. 9: fixed HSD'); save_figure_pair(fig,outbase); close(fig);
end
function value=local_opt(s,name,default), if isfield(s,name), value=s.(name); else, value=default; end, end
function n=local_frames(base,scale), n=max(1,round(base*scale)); end
