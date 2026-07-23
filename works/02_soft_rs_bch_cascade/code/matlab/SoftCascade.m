classdef SoftCascade
% SoftCascade  MATLAB reproduction of Work 02's RS+BCH/PAM4 study.
% It explicitly reuses BCHOSD from Work 01 for the inner BCH/LLOSD path.
methods(Static)
function ensure_deps()
    here=fileparts(mfilename('fullpath')); work=fileparts(fileparts(here)); works=fileparts(work); repo=fileparts(works);
    addpath(fullfile(repo,'works','01_bch_osd_reproduction','code','matlab'));
end
function rs = rs_create(m,k)
    SoftCascade.ensure_deps(); gf=BCHOSD.gf_init(m); n=gf.n; t=floor((n-k)/2); gp=1;
    for i=1:2*t, gp=BCHOSD.poly_mul(gf,gp,[gf.EXP(i+1),1]); end
    rs=struct('gf',gf,'m',m,'n',n,'k',k,'t',t,'d',n-k+1,'g_poly',gp);
end
function c=rs_encode_systematic(rs,msg)
    assert(numel(msg)==rs.k); [~,rem]=BCHOSD.poly_divmod(rs.gf,[zeros(1,rs.n-rs.k),double(msg(:).')],rs.g_poly);
    parity=[rem,zeros(1,rs.n-rs.k-numel(rem))]; c=[parity,double(msg(:).')];
end
function [c,ok]=rs_bm_decode(rs,r)
    gf=rs.gf; r=double(r(:).'); nz=find(r)-1; S=zeros(1,2*rs.t+1);
    if isempty(nz),c=r;ok=true;return,end
    for ii=1:2*rs.t
        s=0; for j=nz, s=bitxor(s,BCHOSD.gf_mul(gf,r(j+1),gf.EXP(mod(ii*j,gf.n)+1))); end; S(ii+1)=s;
    end
    if ~any(S(2:end)),c=r;ok=true;return,end
    L=0; Lam=1; B=1; b=1; shift=1;
    for ni=1:2*rs.t
        delta=S(ni+1); for j=1:L, if j+1<=numel(Lam)&&Lam(j+1)~=0,delta=bitxor(delta,BCHOSD.gf_mul(gf,Lam(j+1),S(ni-j+1)));end,end
        if delta==0,shift=shift+1;continue,end
        co=BCHOSD.gf_div(gf,delta,b); xm=[zeros(1,shift),B]; z=max(numel(Lam),numel(xm)); T=[Lam,zeros(1,z-numel(Lam))];xm=[xm,zeros(1,z-numel(xm))];
        for j=1:z,T(j)=bitxor(T(j),BCHOSD.gf_mul(gf,co,xm(j)));end
        if 2*L<=ni-1,L=ni-L;B=Lam;b=delta;Lam=T;shift=1;else,Lam=T;shift=shift+1;end
    end
    pos=[]; for p=0:rs.n-1
        v=0;for j=0:numel(Lam)-1,if Lam(j+1)~=0,v=bitxor(v,BCHOSD.gf_mul(gf,Lam(j+1),gf.EXP(mod((gf.n-p)*j,gf.n)+1)));end,end
        if v==0,pos(end+1)=p+1;end %#ok<AGROW>
    end
    if numel(pos)~=L||L>rs.t,c=r;ok=false;return,end
    sx=[0,S(2:end)]; omega=zeros(1,2*rs.t+1);
    for i=1:numel(sx),for j=1:numel(Lam),if i+j-1<=2*rs.t+1&&sx(i)~=0&&Lam(j)~=0,omega(i+j-1)=bitxor(omega(i+j-1),BCHOSD.gf_mul(gf,sx(i),Lam(j)));end,end,end
    der=zeros(1,numel(Lam)-1); for i=2:numel(Lam),if mod(i-1,2)==1,der(i-1)=Lam(i);end,end
    c=r; for pp=pos
        p=pp-1; ov=0;dv=0;for i=0:numel(omega)-1,if omega(i+1)~=0,ov=bitxor(ov,BCHOSD.gf_mul(gf,omega(i+1),gf.EXP(mod(i*(gf.n-p),gf.n)+1)));end,end
        for i=0:numel(der)-1,if der(i+1)~=0,dv=bitxor(dv,BCHOSD.gf_mul(gf,der(i+1),gf.EXP(mod(i*(gf.n-p),gf.n)+1)));end,end
        if dv==0,c=r;ok=false;return,end; ev=BCHOSD.gf_mul(gf,gf.EXP(p+1),BCHOSD.gf_div(gf,ov,dv)); c(pp)=bitxor(c(pp),ev);
    end;ok=true;
end
function [c,ok]=rs_lccbr_decode(rs,r,reliability,eta)
    [~,order]=sort(reliability,'ascend'); lrp=order(1:eta); best=[]; bestScore=-inf;
    for z=0:2^eta-1
        q=r;for b=1:eta,if bitget(z,b),q(lrp(b))=bitxor(q(lrp(b)),1);end,end
        [d,valid]=SoftCascade.rs_bm_decode(rs,q);if ~valid,continue,end
        score=sum(reliability(d==r));if score>bestScore,best=d;bestScore=score;end
    end
    if isempty(best),[c,ok]=SoftCascade.rs_bm_decode(rs,r);else,c=best;ok=true;end
end
function x=bits_to_pam4(bits)
    bits=double(bits(:).'); assert(mod(numel(bits),2)==0); x=zeros(1,numel(bits)/2);
    for i=1:numel(x),b1=bits(2*i-1);b0=bits(2*i);if b1==0&&b0==0,x(i)=-3;elseif b1==0,x(i)=-1;elseif b0==1,x(i)=1;else,x(i)=3;end,end
end
function bits=pam4_to_bits_hard(y)
    lev=[-3 -1 1 3]; labels=[0 0;0 1;1 1;1 0];bits=zeros(1,2*numel(y));
    for i=1:numel(y),[~,q]=min(abs(y(i)-lev));bits(2*i-1:2*i)=labels(q,:);end
end
function s=sigma_from_ebn0_pam4(db,rate),s=sqrt(5/(4*rate*10^(db/10)));end
function L=pam4_bit_llr(y,sigma)
    lev=[-3 -1 1 3];labs=[0 0;0 1;1 1;1 0];L=zeros(1,2*numel(y));
    for i=1:numel(y),metric=-(y(i)-lev).^2/(2*sigma^2);for b=1:2
        a=SoftCascade.logsumexp(metric(labs(:,b)==0)); z=SoftCascade.logsumexp(metric(labs(:,b)==1));L(2*i-2+b)=a-z;
    end,end
end
function v=logsumexp(a),u=max(a);v=u+log(sum(exp(a-u)));end
function L=run_channel(bits,db,rate)
    n=numel(bits); if mod(n,2),bits=[bits,0];end;s=SoftCascade.sigma_from_ebn0_pam4(db,rate);L=SoftCascade.pam4_bit_llr(SoftCascade.bits_to_pam4(bits)+s*randn(1,numel(bits)/2),s);L=L(1:n);
end
function codec=cascade_create(m,krs,tbch,tau,eta)
    SoftCascade.ensure_deps(); rs=SoftCascade.rs_create(m,krs);bch=BCHOSD.bch_create(m,tbch);rsBits=m*rs.n;pad=mod(-rsBits,bch.k);blocks=(rsBits+pad)/bch.k;
    codec=struct('m',m,'k_rs',krs,'t_bch',tbch,'llosd_tau',tau,'lcc_eta',eta,'rs',rs,'bch',bch,'rs_bits',rsBits,'n_pad_bits',pad,'n_bch_blocks',blocks,'n_coded_bits',blocks*bch.n,'effective_rate',(krs*m)/(blocks*bch.n));
    [codec.pivots,codec.msgInv]=SoftCascade.bch_recovery(bch);
end
function bits=cascade_encode(codec,msg)
    cr=SoftCascade.rs_encode_systematic(codec.rs,msg);raw=SoftCascade.symbols_to_bits(cr,codec.m);raw=[raw,zeros(1,codec.n_pad_bits)];bits=zeros(1,codec.n_coded_bits);
    for b=1:codec.n_bch_blocks,ii=(b-1)*codec.bch.k+1:b*codec.bch.k;jj=(b-1)*codec.bch.n+1:b*codec.bch.n;bits(jj)=BCHOSD.bch_encode(codec.bch,raw(ii));end
end
function [msg,stats]=cascade_decode(codec,L,scheme)
    if nargin<3,scheme='scheme_a';end; rec=zeros(1,codec.n_bch_blocks*codec.bch.k);rel=zeros(1,codec.n_bch_blocks);
    for b=1:codec.n_bch_blocks,jj=(b-1)*codec.bch.n+1:b*codec.bch.n;[cw,~]=BCHOSD.llosd_decode(codec.bch,L(jj),codec.llosd_tau,true);rec((b-1)*codec.bch.k+1:b*codec.bch.k)=mod(cw(codec.pivots)*codec.msgInv,2);rel(b)=sum(abs(L(jj)));end
    rr=SoftCascade.bits_to_symbols(rec(1:codec.rs_bits),codec.m);sr=zeros(1,codec.rs.n);for i=1:codec.rs.n,ids=floor(((i-1)*codec.m:(i*codec.m-1))/codec.bch.k)+1;sr(i)=mean(rel(ids));end
    [cd,ok]=SoftCascade.rs_lccbr_decode(codec.rs,rr,sr,codec.lcc_eta);if ok,msg=cd(codec.rs.n-codec.rs.k+1:end);else,msg=rr(codec.rs.n-codec.rs.k+1:end);end
    stats=struct('ok',ok,'scheme',scheme);
end
function [p,inv]=bch_recovery(bch)
    G=mod(bch.G,2);[k,n]=size(G);p=[];A=G;
    for col=1:n
        if numel(p)==k,break,end;q=find(A(numel(p)+1:k,col),1);if isempty(q),continue,end;q=q+numel(p);p(end+1)=col;A([numel(p) q],:)=A([q numel(p)],:);for r=1:k,if r~=numel(p)&&A(r,col),A(r,:)=bitxor(A(r,:),A(numel(p),:));end,end
    end
    M=G(:,p);inv=eye(k);for i=1:k,q=find(M(i:k,i),1)+i-1;M([i q],:)=M([q i],:);inv([i q],:)=inv([q i],:);for r=1:k,if r~=i&&M(r,i),M(r,:)=bitxor(M(r,:),M(i,:));inv(r,:)=bitxor(inv(r,:),inv(i,:));end,end,end
end
function b=symbols_to_bits(x,m),b=zeros(1,numel(x)*m);for i=1:numel(x),for j=1:m,b((i-1)*m+j)=bitget(x(i),j);end,end,end
function x=bits_to_symbols(b,m),x=zeros(1,numel(b)/m);for i=1:numel(x),for j=1:m,x(i)=x(i)+b((i-1)*m+j)*2^(j-1);end,end,end
function report=selftest()
    SoftCascade.ensure_deps();rs=SoftCascade.rs_create(5,27);m=mod(0:26,32);cw=SoftCascade.rs_encode_systematic(rs,m);[d,ok]=SoftCascade.rs_bm_decode(rs,cw);codec=SoftCascade.cascade_create(5,27,1,1,2);bits=SoftCascade.cascade_encode(codec,m);[h,st]=SoftCascade.cascade_decode(codec,(1-2*bits)*30,'scheme_a');
    report=struct('name',{'rs_noiseless','cascade_noiseless','pam4_hard'},'pass',{ok&&isequal(d,cw),st.ok&&isequal(h,m),isequal(SoftCascade.pam4_to_bits_hard(SoftCascade.bits_to_pam4([0 0 0 1 1 1 1 0])),[0 0 0 1 1 1 1 0])});
end
end
end
