function out = fig08_rate_sweep(options)
%FIG08_RATE_SWEEP Reduced Fig. 8 rate sweep, retaining the Python defaults.

if nargin<1, options=struct(); end
paths=work01_setup(); visible=local_opt(options,'visible','off');
maxFrames=local_opt(options,'max_frames',800); minErrors=local_opt(options,'min_errors',20);
candidates=[1 4;2 6;3 6;4 8;5 8;7 10]; rates=zeros(1,size(candidates,1));
snrH=zeros(1,6); snrY=zeros(1,6); snrC=zeros(1,6); opsH=zeros(1,6); opsY=zeros(1,6); opsC=zeros(1,6);
for i=1:size(candidates,1)
    t=candidates(i,1); eta=candidates(i,2); code=work01.Core.bch_code(7,t); rates(i)=code.k/code.n;
    fprintf('\n(127,%d), rate=%.3f, t=%d\n',code.k,rates(i),t);
    [snrH(i),opsH(i)]=work01.Experiment.find_snr_for_fer(code,@(c,L)work01.Core.hsd_fast(c,L,1,eta),1e-2,[2,7],maxFrames,minErrors);
    [snrY(i),opsY(i)]=work01.Experiment.find_snr_for_fer(code,@(c,L)work01.Core.ysvl_osd_decode(c,L,1),1e-2,[2,7],maxFrames,minErrors);
    [snrC(i),opsC(i)]=work01.Experiment.find_snr_for_fer(code,@(c,L)work01.Core.cj_osd_decode(c,L,1),1e-2,[2,7],maxFrames,minErrors);
    if isnan(opsH(i)),opsH(i)=1;end; if isnan(opsY(i)),opsY(i)=1;end; if isnan(opsC(i)),opsC(i)=1;end
    fprintf('  HSD(1,%d): %.1f dB, %.0f ops\n',eta,snrH(i),opsH(i));
end
fig=figure('Visible',visible,'Position',[100,100,1100,450]);
subplot(1,2,1); hold on; plot(rates,snrH,'^-k'); plot(rates,snrY,'s-','Color',[.85 .65 0]); plot(rates,snrC,'x-','Color',[1 .5 0]); xlabel('Code rate'); ylabel('Eb/N0 (dB)'); title('(a) Required SNR for FER=10^{-2}'); grid on; legend('HSD','YSVL','CJ');
subplot(1,2,2); hold on; semilogy(rates,opsH,'^-k'); semilogy(rates,opsY,'s-','Color',[.85 .65 0]); semilogy(rates,opsC,'x-','Color',[1 .5 0]); xlabel('Code rate'); ylabel('F_2 / F_{128} ops per frame'); title('(b) Ops for FER=10^{-2}'); grid on; legend('HSD','YSVL','CJ');
save_figure_pair(fig,fullfile(paths.figures_dir,'fig08_rate_sweep')); close(fig);
snr=containers.Map({'HSD','YSVL','CJ'},{snrH,snrY,snrC}); ops=containers.Map({'HSD','YSVL','CJ'},{opsH,opsY,opsC});
out=containers.Map({'rates','snr_1e2','ops_1e2'},{rates,snr,ops}); write_json_file(fullfile(paths.data_dir,'fig08_rate_sweep.json'),out);
end
function value=local_opt(s,name,default), if isfield(s,name), value=s.(name); else, value=default; end, end
