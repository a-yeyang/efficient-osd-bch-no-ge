classdef Work01Experiments
% Work01Experiments  MATLAB counterparts of the Work 01 Python experiments.
% profile=true is a deterministic smoke configuration; profile=false retains
% the published code parameters and is intentionally expensive.
methods(Static)
function out=fig02_nbch(profile)
    if nargin<1,profile=true;end; code=BCHOSD.bch_create(5,2); snr=2:2:8;if profile,ntrials=8;else,ntrials=10000;end;avg=zeros(size(snr));
    for a=1:numel(snr),v=zeros(1,ntrials);for q=1:ntrials,L=(1-2*zeros(1,code.n))+BCHOSD.sigma_from_ebn0(snr(a),code.k/code.n)*randn(1,code.n);[~,s]=BCHOSD.llosd_decode(code,2*L/BCHOSD.sigma_from_ebn0(snr(a),code.k/code.n)^2,2,false);v(q)=s.n_bch_candidates;end;avg(a)=mean(v);end
    out=struct('ebn0_db',snr,'nbch',avg);
end
function out=fer_curve(code,decoder,snr,minErrors,maxFrames)
    out=struct('ebn0_db',snr,'fer',zeros(size(snr)));x=ones(1,code.n);
    for a=1:numel(snr),sigma=BCHOSD.sigma_from_ebn0(snr(a),code.k/code.n);ne=0;nf=0;while nf<maxFrames&&ne<minErrors,L=2*(x+sigma*randn(size(x)))/sigma^2;d=decoder(L);ne=ne+any(d~=0);nf=nf+1;end;out.fer(a)=ne/nf;end
end
function out=fig03_04_fer(profile)
    if nargin<1,profile=true;end;if profile,me=3;mf=20;else,me=60;mf=15000;end;c31=BCHOSD.bch_create(5,2);c63=BCHOSD.bch_create(6,3);out=struct();out.fig03.llosd1=Work01Experiments.fer_curve(c31,@(L) BCHOSD.llosd_decode(c31,L,1,true),2:2:6,me,mf);out.fig04.llosd3=Work01Experiments.fer_curve(c63,@(L) BCHOSD.llosd_decode(c63,L,3,true),3:2:7,me,mf);
end
function out=fig05_comparison(profile)
    if nargin<1,profile=true;end;c=BCHOSD.bch_create(6,3);if profile,me=3;mf=15;else,me=60;mf=15000;end;snr=3:2:7;out=struct('llosd3',Work01Experiments.fer_curve(c,@(L) BCHOSD.llosd_decode(c,L,3,true),snr,me,mf),'sllosd32',Work01Experiments.fer_curve(c,@(L) BCHOSD.sllosd_decode(c,L,[3 2],true),snr,me,mf));
end
function out=fig07_09_longcode(profile)
    if nargin<1,profile=true;end;if profile,me=2;mf=8;else,me=30;mf=3000;end;c127=BCHOSD.bch_create(7,2);c255=BCHOSD.bch_create(8,2);out=struct('fig07',Work01Experiments.fer_curve(c127,@(L) BCHOSD.hsd_decode(c127,L,1,4,true),4:2:8,me,mf),'fig09',Work01Experiments.fer_curve(c255,@(L) BCHOSD.hsd_decode(c255,L,1,4,true),5:2:9,me,mf));
end
function out=fig08_rate_sweep(profile)
    if nargin<1,profile=true;end;ts=[1 2 3];out=struct('t',ts,'hsd_teps',zeros(size(ts)));
    for i=1:numel(ts),c=BCHOSD.bch_create(7,ts(i));[~,s]=BCHOSD.hsd_decode(c,ones(1,c.n)*8,1,4,true);out.hsd_teps(i)=s.n_teps;end
end
function out=fig10_11_ops(profile)
    if nargin<1,profile=true;end;c=BCHOSD.bch_create(5,2);snr=2:2:8;nt=zeros(size(snr));nv=zeros(size(snr));
    for i=1:numel(snr),L=ones(1,c.n)*snr(i);[~,a]=BCHOSD.llosd_decode(c,L,2,true);[~,b]=BCHOSD.hsd_decode(c,L,1,3,true);nt(i)=a.n_teps;nv(i)=b.n_tvs;end;out=struct('ebn0_db',snr,'nteps',nt,'ntvs',nv);
end
function out=tables(profile)
    if nargin<1,profile=true;end;c=BCHOSD.bch_create(5,2);L=ones(1,c.n)*10;[~,o]=BCHOSD.osd_decode(c,L,1,true);[~,l]=BCHOSD.llosd_decode(c,L,2,true);[~,s]=BCHOSD.sllosd_decode(c,L,[2 1],true);[~,h]=BCHOSD.hsd_decode(c,L,1,3,true);out=struct('osd',o,'llosd',l,'sllosd',s,'hsd',h);
end
end
end
