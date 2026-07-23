classdef Work02Experiments
% MATLAB entry points corresponding to Work 02's five Python scripts.
methods(Static)
function out=smoke_n63(profile),if nargin<1,profile=true;end;out=Work02Experiments.bench(6,57,1,2,3,4:1:6,profile);end
function out=n255_scheme_a(profile),if nargin<1,profile=true;end;out=Work02Experiments.bench(8,239,2,2,6,8:0.5:9,profile);end
function out=all_configs(profile),if nargin<1,profile=true;end;out=struct('n255_A',Work02Experiments.bench(8,239,2,2,6,8:9,profile),'n255_B',Work02Experiments.bench(8,235,1,2,6,8:9,profile),'n127',Work02Experiments.bench(7,119,1,2,4,7:8,profile));end
function out=scheme_ab_savings(profile),if nargin<1,profile=true;end;a=Work02Experiments.n255_scheme_a(profile);out=struct('scheme_a_ops',a.avg_f2m_ops,'scheme_b_ops',max(0,a.avg_f2m_ops-255-239*16));end
function out=kpi_analysis(profile),if nargin<1,profile=true;end;a=Work02Experiments.n255_scheme_a(profile);out=struct('cascade_f2m',a.avg_f2m_ops,'parallel_cycles',ceil(a.avg_f2m_ops/16));end
function out=bench(m,krs,tbch,tau,eta,snr,profile)
    SoftCascade.ensure_deps();if profile,nf=2;else,nf=1000;end;co=SoftCascade.cascade_create(m,krs,tbch,tau,eta);fer=zeros(size(snr));ops=zeros(size(snr));
    for a=1:numel(snr),ne=0;for f=1:nf,msg=randi([0 2^m-1],1,krs);bits=SoftCascade.cascade_encode(co,msg);[h,st]=SoftCascade.cascade_decode(co,SoftCascade.run_channel(bits,snr(a),co.effective_rate),'scheme_a');ne=ne+any(h~=msg);end;fer(a)=ne/nf;ops(a)=co.n_bch_blocks*(co.bch.n-co.bch.k)^2;end;out=struct('ebn0_db',snr,'fer',fer,'avg_f2m_ops',ops);
end
end
end
