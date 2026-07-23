classdef BCHOSD
% BCHOSD  MATLAB reproduction of Work 01's finite-field, BCH and OSD code.
% All vectors use MATLAB's 1-based positions; field elements retain the
% Python representation 0..2^m-1.  Polynomial coefficient vectors are low
% degree first, matching src/gf.py.
methods(Static)
function gf = gf_init(m)
    prims = [0 3 7 11 19 37 67 137 285];
    assert(m >= 1 && m <= 8, 'Unsupported GF exponent');
    n = 2^m - 1; expT = zeros(1, 2*n+2); logT = -ones(1,2^m);
    x = 1;
    for i = 0:n-1
        expT(i+1) = x; logT(x+1) = i;
        x = bitshift(x,1);
        if bitand(x, bitshift(1,m)), x = bitxor(x,prims(m+1)); end
    end
    for i=n:2*n+1, expT(i+1)=expT(i-n+1); end
    gf = struct('m',m,'n',n,'prim',prims(m+1),'EXP',expT,'LOG',logT);
end
function c = gf_add(a,b), c = bitxor(a,b); end
function c = gf_mul(gf,a,b)
    a=double(a); b=double(b); c=zeros(size(a+b));
    [a,b] = BCHOSD.broadcast(a,b); nz=(a~=0)&(b~=0);
    c(nz)=gf.EXP(gf.LOG(a(nz)+1)+gf.LOG(b(nz)+1)+1);
end
function c = gf_div(gf,a,b)
    a=double(a); b=double(b); assert(~any(b(:)==0),'division by zero');
    c=zeros(size(a+b)); [a,b]=BCHOSD.broadcast(a,b); nz=a~=0;
    c(nz)=gf.EXP(gf.LOG(a(nz)+1)-gf.LOG(b(nz)+1)+gf.n+1);
end
function c = gf_inv(gf,a)
    assert(~any(a(:)==0),'inverse of zero'); c=gf.EXP(gf.n-gf.LOG(a+1)+1);
end
function c = gf_pow(gf,a,e)
    if a==0, c=double(e==0); return; end
    c=gf.EXP(mod(gf.LOG(a+1)*e,gf.n)+1);
end
function y = poly_eval(gf,p,x)
    y=0; for i=numel(p):-1:1, y=bitxor(BCHOSD.gf_mul(gf,y,x),p(i)); end
end
function out = poly_mul(gf,a,b)
    out=zeros(1,numel(a)+numel(b)-1);
    for i=1:numel(a), for j=1:numel(b)
        if a(i)~=0 && b(j)~=0, out(i+j-1)=bitxor(out(i+j-1),BCHOSD.gf_mul(gf,a(i),b(j))); end
    end, end
end
function [q,r] = poly_divmod(gf,num,den)
    den=den(:).'; while numel(den)>1 && den(end)==0, den(end)=[]; end
    assert(any(den),'zero divisor'); r=num(:).'; q=zeros(1,max(0,numel(r)-numel(den)+1));
    il=BCHOSD.gf_inv(gf,den(end));
    for i=numel(q):-1:1
        if r(i+numel(den)-1)~=0
            co=BCHOSD.gf_mul(gf,r(i+numel(den)-1),il); q(i)=co;
            for j=1:numel(den), r(i+j-1)=bitxor(r(i+j-1),BCHOSD.gf_mul(gf,co,den(j))); end
        end
    end
    while ~isempty(r) && r(end)==0, r(end)=[]; end
end
function g = bch_generator_poly(gf,t)
    seen=[]; g=1;
    for i=1:2*t
        coset=[]; j=mod(i,gf.n);
        while ~ismember(j,coset), coset(end+1)=j; j=mod(2*j,gf.n); end %#ok<AGROW>
        rep=min(coset); if ismember(rep,seen), continue; end; seen(end+1)=rep;
        mp=1; for s=coset, mp=BCHOSD.poly_mul(gf,mp,[gf.EXP(s+1),1]); end
        assert(all(mp==0|mp==1),'minimal polynomial must be binary');
        g=mod(BCHOSD.poly_mul(gf,g,mp),2);
    end
end
function code = bch_create(m,t)
    gf=BCHOSD.gf_init(m); n=gf.n; gp=BCHOSD.bch_generator_poly(gf,t); k=n-numel(gp)+1;
    G=zeros(k,n); for i=1:k, G(i,i:i+numel(gp)-1)=gp; end
    Hext=zeros(m*2*t,n);
    for j=0:n-1, for i=1:2*t
        e=BCHOSD.gf_pow(gf,2,mod(i*j,n));
        for b=0:m-1, Hext((i-1)*m+b+1,j+1)=bitget(e,b+1); end
    end,end
    H=BCHOSD.rref2(Hext); H=H(any(H,2),:);
    assert(all(mod(G*H.',2)==0,'all'),'G H'' mismatch');
    code=struct('m',m,'t',t,'gf',gf,'n',n,'k',k,'g_poly',gp,'G',G,'H',H,'d_design',2*t+1);
end
function c = bch_encode(code,msg), c=mod(double(msg(:).')*code.G,2); end
function [c,ok] = bch_bm_decode(code,r)
    r=double(r(:).'); gf=code.gf; S=zeros(1,2*code.t+1); inds=find(r)-1;
    for i=1:2*code.t
        s=0; for j=inds, s=bitxor(s,BCHOSD.gf_pow(gf,2,mod(i*j,gf.n))); end; S(i+1)=s;
    end
    if ~any(S(2:end)), c=r; ok=true; return; end
    L=0; Lam=1; B=1; b=1; shift=1;
    for ni=1:2*code.t
        delta=S(ni+1); for i=1:L
            if i+1<=numel(Lam) && Lam(i+1)~=0, delta=bitxor(delta,BCHOSD.gf_mul(gf,Lam(i+1),S(ni-i+1))); end
        end
        if delta==0, shift=shift+1; continue; end
        co=BCHOSD.gf_div(gf,delta,b); xmB=[zeros(1,shift),B]; z=max(numel(Lam),numel(xmB));
        T=[Lam,zeros(1,z-numel(Lam))]; xmB=[xmB,zeros(1,z-numel(xmB))];
        for i=1:z, T(i)=bitxor(T(i),BCHOSD.gf_mul(gf,co,xmB(i))); end
        if 2*L<=ni-1, L=ni-L; B=Lam; b=delta; Lam=T; shift=1; else, Lam=T; shift=shift+1; end
    end
    pos=[]; for ii=0:code.n-1
        v=0; for j=0:numel(Lam)-1
            if Lam(j+1)~=0, v=bitxor(v,BCHOSD.gf_mul(gf,Lam(j+1),BCHOSD.gf_pow(gf,2,mod((gf.n-ii)*j,gf.n)))); end
        end
        if v==0, pos(end+1)=ii+1; end %#ok<AGROW>
    end
    if numel(pos)~=L || L>code.t, c=r; ok=false; return; end
    c=r; c(pos)=1-c(pos); ok=true;
end
function [Gsys,colPerm,ops] = gaussian_elim_binary(G)
    G=mod(G,2); [k,n]=size(G); colPerm=1:n; ops=0;
    for i=1:k
        pc=find(G(i,i:n),1); if isempty(pc)
            pr=[]; for rr=i+1:k, if any(G(rr,i:n)), pr=rr; break; end, end
            assert(~isempty(pr),'degenerate generator'); G([i pr],:)=G([pr i],:); pc=find(G(i,i:n),1);
        end
        pc=pc+i-1; if pc~=i, G(:,[i pc])=G(:,[pc i]); colPerm([i pc])=colPerm([pc i]); end
        for rr=1:k, if rr~=i && G(rr,i), G(rr,:)=bitxor(G(rr,:),G(i,:)); ops=ops+n; end, end
    end
    Gsys=G;
end
function [c,stats] = osd_decode(code,L,tau,early)
    if nargin<4, early=true; end; t0=tic; n=code.n; k=code.k; hard=double(L(:).'<0);
    [~,perm]=sort(abs(L),'descend'); [Gs,cp,ops]=BCHOSD.gaussian_elim_binary(code.G(:,perm)); pe=perm(cp); rs=hard(pe); ls=L(pe); bestD=inf; best=[]; count=0;
    for w=0:tau
        sup=BCHOSD.combs(k,w); for q=1:size(sup,1)
            f=rs(1:k); f(sup(q,sup(q,:)~=0))=1-f(sup(q,sup(q,:)~=0)); cand=mod(f*Gs,2); D=BCHOSD.correlation_distance(ls,rs,cand); count=count+1;
            if D<bestD, bestD=D; best=cand; if early && BCHOSD.ml_ok(abs(ls(rs==cand)),code.d_design,sum(rs~=cand),D), break; end, end
        end
        if ~isempty(best) && early && BCHOSD.ml_ok(abs(ls(rs==best)),code.d_design,sum(rs~=best),bestD), break; end
    end
    c=zeros(1,n); c(pe)=best; stats=struct('n_teps',count,'n_bch_candidates',count,'f2',ops,'latency_us',toc(t0)*1e6);
end
function [G,thetaC] = build_rs_systematic_generator(gf,theta,n)
    kp=numel(theta); theta=theta(:).'; thetaC=setdiff(1:n,theta,'stable'); G=zeros(kp,n);
    for i=1:kp, G(i,theta(i))=1; end
    for i=1:kp
        ai=gf.EXP(theta(i)); den=1; for j=1:kp, if j~=i, den=BCHOSD.gf_mul(gf,den,bitxor(ai,gf.EXP(theta(j)))); end,end
        for jj=thetaC
            num=1; ac=gf.EXP(jj); for p=1:kp, if p~=i, num=BCHOSD.gf_mul(gf,num,bitxor(ac,gf.EXP(theta(p)))); end,end
            G(i,jj)=BCHOSD.gf_div(gf,num,den);
        end
    end
end
function [c,stats] = llosd_decode(code,L,tau,early)
    if nargin<4, early=true; end; t0=tic; n=code.n; kp=n-2*code.t; hard=double(L(:).'<0); [~,perm]=sort(abs(L),'descend'); theta=perm(1:kp);
    [G,thetaC]=BCHOSD.build_rs_systematic_generator(code.gf,theta,n); Gpc=G(:,thetaC); u0=hard(theta); vp=zeros(1,numel(thetaC)); for i=find(u0), vp=bitxor(vp,Gpc(i,:)); end
    best=[]; bestD=inf; nt=0; nb=0; terminated=false;
    for w=0:tau
        sup=BCHOSD.combs(kp,w); for q=1:size(sup,1)
            s=sup(q,sup(q,:)~=0); v=vp; for ii=s, v=bitxor(v,Gpc(ii,:)); end; nt=nt+1;
            if any(v>1), continue; end; nb=nb+1; cand=zeros(1,n); cand(theta)=u0; cand(theta(s))=1-cand(theta(s)); cand(thetaC)=v; D=BCHOSD.correlation_distance(L,hard,cand);
            if D<bestD, bestD=D; best=cand; if early && BCHOSD.ml_ok(abs(L(hard==cand)),code.d_design,sum(hard~=cand),D), terminated=true; break; end, end
        end
        if terminated, break; end
    end
    if isempty(best), best=hard; end
    c=best;
    stats=struct('n_teps',nt,'n_bch_candidates',nb,'terminated_early',terminated,'f2m',2*(n*n-kp*kp+kp),'latency_us',toc(t0)*1e6);
end
function [c,stats] = sllosd_decode(code,L,thetaTuple,early)
    if nargin<4, early=true; end; t0=tic; n=code.n; k=code.k; kp=n-2*code.t; hard=double(L(:).'<0); [~,perm]=sort(abs(L),'descend'); theta=perm(1:kp);
    [G,thetaC]=BCHOSD.build_rs_systematic_generator(code.gf,theta,n); u0=hard(theta); vp=zeros(1,numel(thetaC)); for i=find(u0), vp=bitxor(vp,G(i,thetaC)); end
    best=[]; bestD=inf; nt=0; nb=0; terminated=false;
    for rho=0:numel(thetaTuple)-1
        sy=BCHOSD.combs(k,rho); swMax=thetaTuple(rho+1);
        for iy=1:size(sy,1), for rp=0:swMax
            sw=BCHOSD.combs(kp-k,rp)+k;
            for iw=1:size(sw,1)
                s=[sy(iy,sy(iy,:)~=0),sw(iw,sw(iw,:)~=0)]; v=vp; for ii=s, v=bitxor(v,G(ii,thetaC)); end; nt=nt+1;
                if any(v>1), continue; end; nb=nb+1; cand=zeros(1,n); cand(theta)=u0; cand(theta(s))=1-cand(theta(s)); cand(thetaC)=v; D=BCHOSD.correlation_distance(L,hard,cand);
                if D<bestD, bestD=D; best=cand; if early && BCHOSD.ml_ok(abs(L(hard==cand)),code.d_design,sum(hard~=cand),D), terminated=true; break; end,end
            end
            if terminated, break; end
        end,end
        if terminated, break; end
    end
    if isempty(best), best=hard; end; c=best; stats=struct('n_teps',nt,'n_bch_candidates',nb,'latency_us',toc(t0)*1e6);
end
function [c,stats] = hsd_decode(code,L,tau,eta,early)
    if nargin<5, early=true; end; [best,s]=BCHOSD.llosd_decode(code,L,tau,false); hard=double(L(:).'<0); bestD=BCHOSD.correlation_distance(L,hard,best); [~,p]=sort(abs(L),'ascend'); psi=p(1:eta); ntv=0; skip=0;
    for z=0:2^eta-1
        e=zeros(1,code.n); for b=1:eta, e(psi(b))=bitget(z,b); end; r=bitxor(hard,e); ntv=ntv+1;
        if sum(best~=r)<=code.t, skip=skip+1; continue; end
        [d,ok]=BCHOSD.bch_bm_decode(code,r); if ~ok, continue; end; D=BCHOSD.correlation_distance(L,hard,d);
        if D<bestD, best=d; bestD=D; if early && BCHOSD.ml_ok(abs(L(hard==d)),code.d_design,sum(hard~=d),D), break; end,end
    end
    c=best; stats=s; stats.n_tvs=ntv; stats.n_tvs_skipped=skip;
end
function [c,stats] = ysvl_osd_decode(code,L,tau), [c,stats]=BCHOSD.osd_decode(code,L,tau,true); end
function [c,stats] = cj_osd_decode(code,L,tau), [c,stats]=BCHOSD.osd_decode(code,L,tau,true); end
function [c,stats] = plcc_decode(code,L,eta), [c,stats]=BCHOSD.hsd_decode(code,L,0,eta,true); end
function c = ml_decode_full_codebook(code,L)
    assert(code.k<=20,'Full ML only intended for short codes'); bestD=inf; c=[];
    for z=0:2^code.k-1
        m=zeros(1,code.k); for b=1:code.k,m(b)=bitget(z,b);end; x=BCHOSD.bch_encode(code,m); D=BCHOSD.correlation_distance(L,double(L(:).'<0),x); if D<bestD,bestD=D;c=x;end
    end
end
function r = run_mc(code,decoder,ebn0,nFrames,seed)
    rng(seed); r=struct('ebn0_db',ebn0,'fer',zeros(size(ebn0)));
    zero=zeros(1,code.n); x=1-2*zero;
    for a=1:numel(ebn0), sig=BCHOSD.sigma_from_ebn0(ebn0(a),code.k/code.n); ne=0;
        for f=1:nFrames, L=2*(x+sig*randn(size(x)))/(sig^2); d=decoder(L); ne=ne+any(d~=zero); end
        r.fer(a)=ne/nFrames;
    end
end
function x = bpsk_modulate(c), x=1-2*double(c); end
function y = awgn(x,sigma), y=x+sigma*randn(size(x)); end
function s = sigma_from_ebn0(db,rate), s=sqrt(1/(2*rate*10^(db/10))); end
function L = llr_from_y(y,s), L=2*y/(s*s); end
function D = correlation_distance(L,r,c), D=sum(abs(L(r~=c))); end
function ok = ml_ok(matchAbs,d,dw,D)
    K=max(0,d-dw-1); matchAbs=sort(matchAbs); ok=(K==0 && D<=0) || (K>0 && D<=sum(matchAbs(1:min(K,numel(matchAbs)))));
end
function c = combs(n,w)
    if w==0, c=zeros(1,0); else, c=nchoosek(1:n,w); end
end
function A = rref2(A)
    A=mod(A,2); [r,cols]=size(A); p=1; for j=1:cols
        q=find(A(p:r,j),1); if isempty(q),continue,end; q=q+p-1; A([p q],:)=A([q p],:);
        for i=1:r, if i~=p && A(i,j), A(i,:)=bitxor(A(i,:),A(p,:)); end,end; p=p+1; if p>r,break,end
    end
end
function [a,b]=broadcast(a,b)
    if isscalar(a),a=repmat(a,size(b)); elseif isscalar(b),b=repmat(b,size(a)); else, assert(isequal(size(a),size(b))); end
end
function report = selftest()
    code=BCHOSD.bch_create(5,2); msg=mod(0:code.k-1,2); cw=BCHOSD.bch_encode(code,msg); [bm,ok]=BCHOSD.bch_bm_decode(code,cw); L=(1-2*cw)*12;
    [o,~]=BCHOSD.osd_decode(code,L,1,true); [l,~]=BCHOSD.llosd_decode(code,L,1,true); [s,~]=BCHOSD.sllosd_decode(code,L,[2 1],true); [h,~]=BCHOSD.hsd_decode(code,L,1,3,true);
    report=struct('name',{'encode_bm','osd','llosd','sllosd','hsd'},'pass',{ok&&isequal(bm,cw),isequal(o,cw),isequal(l,cw),isequal(s,cw),isequal(h,cw)});
end
end
end
