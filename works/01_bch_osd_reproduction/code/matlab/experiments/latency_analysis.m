function latency_analysis(options)
%LATENCY_ANALYSIS Plot latency figures from MATLAB-generated Tables I/III/IV.

if nargin<1, options=struct(); end
paths=work01_setup(); visible=local_opt(options,'visible','off');
t1=jsondecode(fileread(fullfile(paths.tables_dir,'table_I.json')));
t3=jsondecode(fileread(fullfile(paths.tables_dir,'table_III.json')));
t4=jsondecode(fileread(fullfile(paths.tables_dir,'table_IV.json'))); %#ok<NASGU>
snrKeys={'4dB','5dB','6dB'}; snrVals=[4,5,6];
data1=local_extract(t1,{'OSD(1)','LLOSD(3)','LLOSD-B(3)'},snrKeys);
data2=local_extract(t3,{'OSD(1)','LLOSD(3)','HSD(1,6)'},snrKeys);
fig=figure('Visible',visible,'Position',[100,100,1100,450]);
subplot(1,2,1); bar(data1.'); set(gca,'YScale','log','XTickLabel',snrKeys); ylabel('Latency (us)'); title('(63,45) BCH'); grid on; legend({'OSD(1)','LLOSD(3)','LLOSD-B(3)'});
subplot(1,2,2); bar(data2.'); set(gca,'YScale','log','XTickLabel',snrKeys); ylabel('Latency (us)'); title('(127,99) BCH'); grid on; legend({'OSD(1)','LLOSD(3)','HSD(1,6)'});
save_figure_pair(fig,fullfile(paths.figures_dir,'latency_bars')); close(fig);
fig=figure('Visible',visible,'Position',[100,100,700,450]); hold on;
semilogy(snrVals,data1(1,:),'o-','DisplayName','OSD(1)'); semilogy(snrVals,data1(2,:),'s-','DisplayName','LLOSD(3)'); semilogy(snrVals,data1(3,:),'^-','DisplayName','LLOSD-B(3)');
xlabel('Eb/N0 (dB)'); ylabel('Per-frame decoding latency (us)'); title('(63,45) BCH latency'); grid on; legend;
save_figure_pair(fig,fullfile(paths.figures_dir,'latency_vs_snr')); close(fig);
end

function data=local_extract(decoded,algorithms,snrKeys)
data=zeros(numel(algorithms),numel(snrKeys));
for a=1:numel(algorithms)
    algoField=matlab.lang.makeValidName(algorithms{a}); row=decoded.(algoField);
    for s=1:numel(snrKeys), metric=row.(matlab.lang.makeValidName(snrKeys{s})); data(a,s)=metric.latency_us; end
end
end
function value=local_opt(s,name,default),if isfield(s,name),value=s.(name);else,value=default;end,end
