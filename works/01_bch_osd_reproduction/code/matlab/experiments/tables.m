function allTables = tables(options)
%TABLES Reproduce Tables I-IV, including JSON and Markdown output.
% `n_trials_scale` is provided only for an intentionally short run.

if nargin<1, options=struct(); end
paths=work01_setup(); scale=local_opt(options,'n_trials_scale',1); ebn0s=[4,5,6];

code63=work01.Core.bch_code(6,3);
allTables=struct();
fprintf('=== Table I: (63,45) ===\n');
allTables.table_I=local_build(code63,{ ...
    'OSD(1)',@(c,L)work01.Core.osd_decode(c,L,1); ...
    'LLOSD(3)',@(c,L)work01.Core.llosd_fast(c,L,3,false); ...
    'LLOSD-B(3)',@(c,L)work01.Core.llosd_fast(c,L,3,true)},ebn0s,local_n(500,scale));
local_write_table(paths,'table_I',allTables.table_I,ebn0s,'Table I. Numerical Results in Decoding the (63,45) BCH Code');

fprintf('=== Table II: (63,45) ===\n');
allTables.table_II=local_build(code63,{ ...
    'LLOSD-B(3)',@(c,L)work01.Core.llosd_fast(c,L,3,true); ...
    'SLLOSD-B(3,2)',@(c,L)work01.Core.sllosd_fast(c,L,[3,2],true); ...
    'YSVL OSD(1)',@(c,L)work01.Core.ysvl_osd_decode(c,L,1); ...
    'CJ OSD(1)',@(c,L)work01.Core.cj_osd_decode(c,L,1)},ebn0s,local_n(500,scale));
local_write_table(paths,'table_II',allTables.table_II,ebn0s,'Table II. Numerical Results in Decoding the (63,45) BCH Code');

code127=work01.Core.bch_code(7,4); fprintf('=== Table III: (127,99) ===\n');
allTables.table_III=local_build(code127,{ ...
    'OSD(1)',@(c,L)work01.Core.osd_decode(c,L,1); ...
    'LLOSD(3)',@(c,L)work01.Core.llosd_fast(c,L,3,true); ...
    'HSD(1,8)',@(c,L)work01.Core.hsd_fast(c,L,1,8); ...
    'HSD(1,6)',@(c,L)work01.Core.hsd_fast(c,L,1,6)},ebn0s,local_n(200,scale));
local_write_table(paths,'table_III',allTables.table_III,ebn0s,'Table III. Numerical Results in Decoding the (127,99) BCH Code');

code255=work01.Core.bch_code(8,4); fprintf('=== Table IV: (255,223) ===\n');
allTables.table_IV=local_build(code255,{ ...
    'OSD(1)',@(c,L)work01.Core.osd_decode(c,L,1); ...
    'LLOSD(2)',@(c,L)work01.Core.llosd_fast(c,L,2,true); ...
    'PLCC(8)',@(c,L)work01.Core.plcc_decode(c,L,8); ...
    'HSD(1,8)',@(c,L)work01.Core.hsd_fast(c,L,1,8)},ebn0s,local_n(100,scale));
local_write_table(paths,'table_IV',allTables.table_IV,ebn0s,'Table IV. Numerical Results in Decoding the (255,223) BCH Code');
end

function tbl=local_build(code,pairs,ebn0s,nTrials)
tbl=containers.Map;
for i=1:size(pairs,1)
    name=pairs{i,1}; fn=pairs{i,2}; row=containers.Map;
    for ebn0=ebn0s
        metrics=work01.Experiment.instrument(code,fn,ebn0,nTrials,0);
        row(sprintf('%ddB',ebn0))=metrics;
        fprintf('  %s @ %d dB: F2=%.2e F2m=%.2e FP=%.1f lat=%.1f us\n',name,ebn0,metrics.F2,metrics.F2m,metrics.FP,metrics.latency_us);
    end
    tbl(name)=row;
end
end

function local_write_table(paths,name,tbl,ebn0s,titleText)
write_json_file(fullfile(paths.tables_dir,[name '.json']),tbl);
fid=fopen(fullfile(paths.tables_dir,[name '.md']),'w','n','UTF-8'); cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# %s\n\n',titleText); fprintf(fid,'| Algorithm | Eb/N0 (dB) | F_2 ops | F_2^m ops | Floating ops | Latency (us) |\n'); fprintf(fid,'|---|---|---|---|---|---|\n');
names=tbl.keys;
for i=1:numel(names)
    row=tbl(names{i});
    for j=1:numel(ebn0s)
        metric=row(sprintf('%ddB',ebn0s(j))); algo=''; if j==1,algo=names{i};end
        fprintf(fid,'| %s | %d | %s | %s | %.1f | %.1f |\n',algo,ebn0s(j),local_sci(metric.F2),local_sci(metric.F2m),metric.FP,metric.latency_us);
    end
end
end

function s=local_sci(x), if x<1,s=sprintf('%.2f',x);else,s=sprintf('%.2e',x);end,end
function value=local_opt(s,name,default),if isfield(s,name),value=s.(name);else,value=default;end,end
function n=local_n(base,scale),n=max(1,round(base*scale));end
