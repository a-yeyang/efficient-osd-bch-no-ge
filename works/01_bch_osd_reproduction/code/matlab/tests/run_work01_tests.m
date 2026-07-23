% RUN_WORK01_TESTS Fast deterministic assertions for the Base-MATLAB port.
% Run with: run(fullfile(work01_paths().matlab_root,'tests','run_work01_tests.m'))

testDir = fileparts(mfilename('fullpath'));
matlabRoot = fileparts(testDir);
addpath(genpath(matlabRoot));
paths = work01_setup();
fprintf('Work 01 MATLAB tests: %s\n', matlabRoot);

% Paths must be derived from source location, never the current folder.
previousDir = pwd;
cleanupDir = onCleanup(@() cd(previousDir)); %#ok<NASGU>
cd(tempdir);
derived = work01_paths();
assert(strcmp(derived.work_root, paths.work_root));
assert(isfolder(derived.assets_root));
cd(previousDir);

% GF tables, scalar/vector operations, and polynomial arithmetic.
gf = work01.Core.gf_create(5);
assert(gf.n == 31 && gf.prim == 37);
assert(work01.Core.gf_mul(gf,2,16) == 5);
assert(work01.Core.gf_div(gf,5,16) == 2);
assert(work01.Core.gf_pow(gf,2,5) == 5);
assert(isequal(work01.Core.gf_vmul(gf,[2,0,16],[16,4,0]),[5,0,0]));
assert(isequal(work01.Core.gf_poly_mul(gf,[1,1],[1,1]),[1,0,1]));
[q,r] = work01.Core.gf_poly_divmod(gf,[1,0,1],[1,1]);
assert(isequal(q,[1,1]) && isempty(r));

% Primitive BCH dimensions, parity check, encoder, and BM correction.
code31 = work01.Core.bch_code(5,2);
code63 = work01.Core.bch_code(6,3);
code127 = work01.Core.bch_code(7,4);
code255 = work01.Core.bch_code(8,4);
assert(isequal([code31.n,code31.k],[31,21]));
assert(isequal([code63.n,code63.k],[63,45]));
assert(isequal([code127.n,code127.k],[127,99]));
assert(isequal([code255.n,code255.k],[255,223]));
assert(all(mod(code31.G*code31.H.',2)==0,'all'));
msg = mod(0:code31.k-1,2); cw = work01.Core.bch_encode(code31,msg);
assert(work01.Core.is_codeword(code31,cw));
received = cw; received([1,7]) = bitxor(received([1,7]),1);
[decoded,ok] = work01.Core.bm_decode(code31,received);
assert(ok && isequal(decoded,cw));

% Channel convention and OSD GE systematic form.
assert(isequal(work01.Core.bpsk_modulate([0,1,0]),[1,-1,1]));
assert(abs(work01.Core.sigma_from_ebn0(0,0.5)-1) < 1e-12);
[gSys,~,~] = work01.Core.gaussian_elim_binary(code31.G(:,end:-1:1));
assert(isequal(gSys(:,1:code31.k),eye(code31.k)));

% RS/Lagrange construction plus all decoder families on a known all-zero word.
kPrime = code31.n - 2*code31.t;
[gRS,thetaC] = work01.Core.build_rs_systematic_generator(code31.gf,1:kPrime,kPrime,code31.n);
assert(isequal(gRS(:,1:kPrime),eye(kPrime)) && numel(thetaC)==2*code31.t);
Lgood = 25*ones(1,code31.n);
[o,so] = work01.Core.osd_decode(code31,Lgood,1,true);
[l,sl] = work01.Core.llosd_decode(code31,Lgood,1,false,true);
[lf,slf] = work01.Core.llosd_fast(code31,Lgood,1,true,true);
[seg,ss] = work01.Core.sllosd_fast(code31,Lgood,[1,1],true,true);
[h,sh] = work01.Core.hsd_fast(code31,Lgood,1,2,true,true);
[y,sy] = work01.Core.ysvl_osd_decode(code31,Lgood,1);
[cj,scj] = work01.Core.cj_osd_decode(code31,Lgood,1);
[p,sp] = work01.Core.plcc_decode(code31,Lgood,2);
assert(all(o==0) && all(l==0) && all(lf==0) && all(seg==0) && all(h==0));
assert(all(y==0) && all(cj==0) && all(p==0));
assert(so.n_teps>=1 && sl.n_bch_candidates>=1 && slf.n_teps>=1 && ss.n_teps>=1);
assert(sh.n_teps_llosd>=1 && sp.n_tvs>=1 && sy.counters.fp==so.counters.fp+400);
assert(scj.counters.f2<=so.counters.f2);
assert(work01.Core.sllosd_n_teps_theoretical(45,57,[3,2]) == 3854);

% Cross-language regression vector: obtained from the unmodified Python
% implementation with a no-tie LLR vector.  It exercises nontrivial TEPs,
% BCH filtering, segmented enumeration, Chase/BM, and operation counts.
Lfixed = [-7.66666666667,-1.999,3.66866666667,-6.33033333333,-0.662666666667, ...
    5.005,-4.994,0.673666666667,6.34133333333,-3.65766666667,2.01, ...
    7.67766666667,-2.32133333333,3.34633333333,-6.65266666667,-0.985, ...
    4.68266666667,-5.31633333333,0.351333333333,6.019,-3.98,1.68766666667, ...
    7.35533333333,-2.64366666667,3.024,-6.975,-1.30733333333,4.36033333333, ...
    -5.63866666667,0.029,5.69666666667];
expected = double('1101111001001011011010010110110') - double('0');
[xo,sto] = work01.Core.osd_decode(code31,Lfixed,1,false);
[xl,stl] = work01.Core.llosd_decode(code31,Lfixed,1,false,false);
[xs,sts] = work01.Core.sllosd_fast(code31,Lfixed,[1,1],true,false);
[xh,sth] = work01.Core.hsd_fast(code31,Lfixed,1,2,true,false);
assert(isequal(xo,expected) && sto.n_teps==22);
assert(isequal(xl,expected) && stl.n_teps==28 && stl.n_bch_candidates==1);
assert(isequal(xs,expected) && sts.n_teps==154 && sts.n_bch_candidates==4);
assert(isequal(xh,expected) && sth.n_teps_llosd==28 && sth.n_tvs==4);

% Exact ML on a short code and approximate-ML API smoke test.
code7 = work01.Core.bch_code(3,1);
target = work01.Core.bch_encode(code7,[1,0,1,0]);
[ml,~] = work01.Core.ml_decode_full_codebook(code7,50*work01.Core.bpsk_modulate(target));
assert(isequal(ml,target));
[mla,~] = work01.Core.ml_approx_by_high_order_llosd(code7,50*work01.Core.bpsk_modulate(target),2);
assert(work01.Core.is_codeword(code7,mla));

% Simulator / analytical bound smoke tests use actual decoder calls.
res = work01.Core.run_mc(code7,@(c,L) work01.Core.llosd_fast(c,L,1,true,true), ...
    [3,4],3,99,7,false,0);
assert(numel(res.fer)==2 && all(res.n_frames==3) && all(res.avg_n_teps>=1));
bound = work01.Core.sphere_packing_bound_fer(code31.n,code31.k,[0,3,6]);
assert(all(isfinite(bound)) && bound(1)>=bound(2) && bound(2)>=bound(3));

fprintf('PASS: GF, BCH/BM, OSD/LLOSD/SLLOSD/HSD/baselines/ML/simulator assertions succeeded.\n');
