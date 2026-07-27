function main(fast, useParallel)
%MAIN  One-click RS+BCH low-latency cascade simulation (MATLAB).
%
%   低时延级联码算法仿真 —— MATLAB 实现（对应 Python 版 rs_bch_cascade/）。
%
%   Runs the full pipeline for two configs (n=127, n=255) and four methods:
%     1. Pure RS-BM              (hard-decision baseline)
%     2. Pure RS-LCC-BR          (soft-decision baseline)
%     3. Cascade RS+BCH, inner = LLOSD (Lagrange, THE PAPER'S METHOD) + LCC-BR
%     4. Cascade RS+BCH, inner = OSD   (traditional Gaussian elim, CONTRAST)  + LCC-BR
%
%   Chain per frame: RS+BCH cascade encode -> PAM4 Gray -> AWGN -> per-bit LLR
%   -> inner soft decode -> RS outer soft decode.
%
%   Outputs:
%     figures/matlab_n127_ber.{png,fig}, figures/matlab_n255_ber.{png,fig}
%     figures/matlab_latency_bars.{png,fig}   (op-count KPI, +10% threshold)
%     data/matlab_results.mat
%     Console: per-SNR BER/FER, per-group latency (us) + theoretical F2^m ops,
%              and the "RS+BCH vs pure-RS latency increase <= 10%" KPI verdict.
%
%   PARALLEL: frames are independent Monte-Carlo trials, so the whole sweep is
%   parallelized at the frame level via parfor (Parallel Computing Toolbox).
%   A single Processes parpool is started once and reused across all 8 groups.
%   Each frame draws from its own mrg32k3a substream keyed by the global frame
%   index -> results are reproducible regardless of worker count / scheduling.
%
%   Usage:
%     main               % FAST preset ON, parallel ON if PCT available
%     main(true)         % FAST preset ON
%     main(false)        % full precision (slower)
%     main(true, false)  % FAST but FORCE serial (for speedup comparison)

    if nargin < 1, fast = true; end
    if nargin < 2 || isempty(useParallel)
        useParallel = hasParallel();
    end

    % Make sure this file's folder is on the path (so classes resolve).
    here = fileparts(mfilename('fullpath'));
    addpath(here);
    if ~exist(fullfile(here, 'figures'), 'dir'), mkdir(fullfile(here, 'figures')); end
    if ~exist(fullfile(here, 'data'), 'dir'), mkdir(fullfile(here, 'data')); end

    fprintf('\n############################################################\n');
    fprintf('# RS+BCH 低时延级联码仿真 (MATLAB)\n');
    fprintf('#   创新点: BCH 是 RS 的二元子码 -> 用 Lagrange 插值构造\n');
    fprintf('#   RS 系统生成矩阵 (免高斯消元), 过滤非二元候选 (Theorem 2)\n');
    fprintf('#   -> 降低 OSD 复杂度与时延.\n');
    if fast
        fprintf('#   预设: FAST (数分钟). 关闭用 main(false) 跑完整精度.\n');
    else
        fprintf('#   预设: FULL 精度 (较慢).\n');
    end

    % ---- parallel pool (start once, reuse across all 8 groups) ----
    ncores = feature('numcores');
    if useParallel
        nworkers = start_parpool(ncores);
        if nworkers < 1
            useParallel = false;   % pool failed to start -> fall back to serial
        end
    end
    if useParallel
        fprintf('#   并行: 开 (parfor, %d worker / %d 核; 帧级 mrg32k3a 子流, 可复现)\n', ...
            nworkers, ncores);
    else
        fprintf('#   并行: 关 (串行; %d 核可用)\n', ncores);
    end
    fprintf('############################################################\n\n');

    % ---- numerical self-test first ----
    selftest();

    % ---- configs ----
    % n=127: RS(127,119)+BCH(127,120)  m=7 t=1 tau=1 eta=4
    % n=255: RS(255,239)+BCH(255,239)  m=8 t=2 tau=2 eta=4
    cfgs = struct( ...
        'n127', CascadeConfig(7, 119, 1, 1, 4), ...
        'n255', CascadeConfig(8, 239, 2, 2, 4));

    if fast
        ebn0_127 = 5.0:1.0:10.0;
        ebn0_255 = 6.0:1.0:11.0;
        frames_pure = 120; frames_casc = 60;
    else
        ebn0_127 = 5.0:0.5:11.0;
        ebn0_255 = 6.0:0.5:12.0;
        frames_pure = 400; frames_casc = 200;
    end

    all_results = struct();
    cfg_names = {'n127', 'n255'};
    ebn0_lists = {ebn0_127, ebn0_255};

    t_all = tic;
    for ci = 1:numel(cfg_names)
        name = cfg_names{ci};
        cfg = cfgs.(name);
        ebn0_list = ebn0_lists{ci};
        fprintf('\n######### %s: %s #########\n', name, cfg.describe());

        pure = PureRSCodec(cfg);
        casc = CascadedCodec(cfg);

        methods = {
            'pure_rs_bm',    'Pure RS-BM',                @(msg) pure.encode(msg), @(llr,c) pure.decode(llr,'hard',c),  pure.effective_rate, frames_pure;
            'pure_rs_lccbr', 'Pure RS-LCC-BR',            @(msg) pure.encode(msg), @(llr,c) pure.decode(llr,'soft',c),  pure.effective_rate, frames_pure;
            'cascade_llosd', 'Cascade LLOSD(Lagrange)',   @(msg) casc.encode(msg), @(llr,c) casc.decode(llr,'llosd',c), casc.effective_rate, frames_casc;
            'cascade_osd',   'Cascade OSD(GaussElim)',    @(msg) casc.encode(msg), @(llr,c) casc.decode(llr,'osd',c),   casc.effective_rate, frames_casc;
        };

        results_by_method = struct();
        for mi = 1:size(methods, 1)
            key    = methods{mi, 1};
            label  = methods{mi, 2};
            enc    = methods{mi, 3};
            dec    = methods{mi, 4};
            rate   = methods{mi, 5};
            frames = methods{mi, 6};
            fprintf('\n--- %s / %s ---\n', name, label);
            opts = struct('seed', 0, 'min_frame_errors', 15, ...
                'max_frames', frames, 'verbose', true, 'parallel', useParallel);
            res = run_bench(label, enc, dec, rate, ebn0_list, ...
                cfg.k_rs, cfg.m, opts);
            results_by_method.(key) = res;
        end
        all_results.(name) = results_by_method;

        % ---- BER-SNR semilogy figure ----
        plotBer(name, cfg, results_by_method, here);
    end

    % ---- latency KPI bar chart (uses n=255 config) ----
    plotLatencyKpi(all_results.n255, here);
    printKpiTable(all_results);

    total_wall = toc(t_all);

    save(fullfile(here, 'data', 'matlab_results.mat'), 'all_results');
    fprintf('\n已保存: data/matlab_results.mat, figures/matlab_n127_ber.*, ');
    fprintf('figures/matlab_n255_ber.*, figures/matlab_latency_bars.*\n');
    if useParallel
        fprintf('总墙钟(并行, %d worker): %.1f s   —— 关并行对比: main(%s, false)\n', ...
            nworkers, total_wall, mat2str(fast));
    else
        fprintf('总墙钟(串行): %.1f s   —— 开并行对比: main(%s, true)\n', ...
            total_wall, mat2str(fast));
    end
end


% ======================================================================
function tf = hasParallel()
    % True if the Parallel Computing Toolbox is installed AND licensed.
    tf = ~isempty(ver('parallel')) && license('test', 'Distrib_Computing_Toolbox');
end


% ======================================================================
function nworkers = start_parpool(ncores)
    % Start (or reuse) a Processes parpool. Returns worker count, 0 on failure.
    nworkers = 0;
    try
        pool = gcp('nocreate');
        if isempty(pool)
            % Leave a little headroom; cap at 12 (M4 Pro).
            n = min(12, max(1, ncores));
            pool = parpool('Processes', n);
        end
        nworkers = pool.NumWorkers;
    catch ME
        fprintf('#   [警告] parpool 启动失败, 回退串行: %s\n', ME.message);
        nworkers = 0;
    end
end


% ======================================================================
function plotBer(name, cfg, rbm, here)
    styles = {
        'pure_rs_bm',    'D-',  [0.85 0.1 0.1], 'Pure RS-BM';
        'pure_rs_lccbr', 's-',  [0.9 0.55 0.0], 'Pure RS-LCC-BR';
        'cascade_llosd', '^-',  [0.0 0.3 0.85], 'Cascade LLOSD (Lagrange)';
        'cascade_osd',   'v--', [0.0 0.6 0.3],  'Cascade OSD (Gauss elim)';
    };
    f = figure('Visible', 'off', 'Position', [100 100 760 500]);
    hold on;
    for i = 1:size(styles, 1)
        key = styles{i, 1};
        if ~isfield(rbm, key), continue; end
        r = rbm.(key);
        ber = max(r.ber, 1e-6);
        semilogy(r.ebn0_db, ber, styles{i, 2}, 'Color', styles{i, 3}, ...
            'LineWidth', 1.6, 'MarkerSize', 7, 'DisplayName', styles{i, 4});
    end
    set(gca, 'YScale', 'log');
    grid on; box on;
    xlabel('Eb/N0 (dB)'); ylabel('BER (info bits)');
    title(sprintf('%s: %s', name, cfg.describe()), 'Interpreter', 'none', 'FontSize', 11);
    legend('Location', 'southwest', 'FontSize', 9);
    saveas(f, fullfile(here, 'figures', sprintf('matlab_%s_ber.png', name)));
    savefig(f, fullfile(here, 'figures', sprintf('matlab_%s_ber.fig', name)));
    close(f);
end


% ======================================================================
function plotLatencyKpi(rbm, here)
    % Latency-increase KPI using theoretical F2^m op counts at converged SNR.
    % Case (a): cascade LLOSD vs Pure RS-BM ; Case (b): vs Pure RS-LCC-BR.
    % Also plots cascade OSD to show LLOSD (Lagrange) is cheaper.
    ebn0 = rbm.pure_rs_bm.ebn0_db;
    % Common length across all four methods (cascade methods may early-stop).
    common = min([numel(rbm.pure_rs_bm.avg_f2m_ops), ...
                  numel(rbm.pure_rs_lccbr.avg_f2m_ops), ...
                  numel(rbm.cascade_llosd.avg_f2m_ops), ...
                  numel(rbm.cascade_osd.avg_f2m_ops)]);
    % pick up to 3 high-SNR points where cascade FER has converged (< 0.5)
    casc_fer = rbm.cascade_llosd.fer;
    npts = min(common, numel(casc_fer));
    cand = find(casc_fer(1:npts) < 0.5);
    if isempty(cand), cand = npts; end
    key_idx = cand(max(1, numel(cand)-2):end);

    lcc  = rbm.pure_rs_lccbr.avg_f2m_ops;
    cllo = rbm.cascade_llosd.avg_f2m_ops;
    cosd = rbm.cascade_osd.avg_f2m_ops;

    % Total elementary ops = GF(2^m) ops (f2m) + binary ops (f2).
    % The paper's point: OSD pays a large binary Gaussian-elimination cost
    % (lands in f2), which LLOSD's Lagrange construction avoids (f2 ~ 0).
    % So the fair latency proxy is f2m + f2, NOT f2m alone.
    cllo_tot = cllo + rbm.cascade_llosd.avg_f2_ops;
    cosd_tot = cosd + rbm.cascade_osd.avg_f2_ops;

    snrs = ebn0(key_idx);
    ra = 100 * (cllo(key_idx) - lcc(key_idx)) ./ max(lcc(key_idx), 1);  % vs soft
    ro = 100 * (cosd(key_idx) - lcc(key_idx)) ./ max(lcc(key_idx), 1);

    f = figure('Visible', 'off', 'Position', [100 100 900 420]);
    subplot(1, 2, 1);
    b = bar(categorical(arrayfun(@(s) sprintf('%.1f dB', s), snrs, 'uni', 0)), ...
        [ra(:), ro(:)]);
    b(1).FaceColor = [0.0 0.6 0.3]; b(2).FaceColor = [0.6 0.6 0.6];
    hold on; yline(10, 'k--', '+10% KPI', 'LineWidth', 1.2);
    ylabel('时延增幅 vs 纯 RS-LCC-BR (%)');
    title('级联时延增幅 (F_2^m 运算量)');
    legend({'LLOSD (Lagrange)', 'OSD (Gauss elim)'}, 'Location', 'best', 'FontSize', 8);
    grid on;

    subplot(1, 2, 2);
    b2 = bar(categorical(arrayfun(@(s) sprintf('%.1f dB', s), snrs, 'uni', 0)), ...
        [cllo_tot(key_idx)', cosd_tot(key_idx)']);
    b2(1).FaceColor = [0.0 0.3 0.85]; b2(2).FaceColor = [0.6 0.6 0.6];
    set(gca, 'YScale', 'log');   % OSD (GE) is orders of magnitude larger
    ylabel('平均总运算量 / 帧 (F_2^m + F_2, log)');
    title({'LLOSD(Lagrange) << OSD(高斯消元)', '(对数轴)'}, 'FontSize', 10);
    legend({'LLOSD (Lagrange)', 'OSD (Gauss elim)'}, 'Location', 'northwest', 'FontSize', 8);
    grid on;
    % Annotate the speedup factor above each pair.
    for j = 1:numel(key_idx)
        idx = key_idx(j);
        spd = cosd_tot(idx) / max(cllo_tot(idx), 1);
        text(j, cosd_tot(idx) * 1.3, sprintf('%.0f\\times', spd), ...
            'HorizontalAlignment', 'center', 'FontSize', 9, ...
            'FontWeight', 'bold', 'Color', [0.0 0.3 0.85]);
    end

    saveas(f, fullfile(here, 'figures', 'matlab_latency_bars.png'));
    savefig(f, fullfile(here, 'figures', 'matlab_latency_bars.fig'));
    close(f);
end


% ======================================================================
function printKpiTable(all_results)
    fprintf('\n\n================= 时延 / 运算量 KPI =================\n');
    for name = {'n127', 'n255'}
        nm = name{1};
        if ~isfield(all_results, nm), continue; end
        rbm = all_results.(nm);
        fprintf('\n--- %s ---\n', nm);
        fprintf('%-8s %10s %10s %10s %10s %12s\n', 'Eb/N0', 'RS-BM', ...
            'RS-LCC', 'LLOSD', 'OSD(f2m)', 'OSD(f2)');
        ebn0 = rbm.pure_rs_bm.ebn0_db;
        npts = numel(ebn0);
        for i = 1:npts
            v = @(k) valAt(rbm, k, i, 'avg_f2m_ops');
            osd_f2 = valAt(rbm, 'cascade_osd', i, 'avg_f2_ops');
            fprintf('%-8.2f %10.0f %10.0f %10.0f %10.0f %12.0f\n', ebn0(i), ...
                v('pure_rs_bm'), v('pure_rs_lccbr'), v('cascade_llosd'), ...
                v('cascade_osd'), osd_f2);
        end
        % KPI verdict at the highest SNR that ALL four methods reached
        % (cascade methods may early-terminate at a lower SNR than pure RS).
        lcc  = rbm.pure_rs_lccbr.avg_f2m_ops;
        cllo = rbm.cascade_llosd.avg_f2m_ops;
        cllo_f2 = rbm.cascade_llosd.avg_f2_ops;
        cosd_f2m = rbm.cascade_osd.avg_f2m_ops;
        cosd_f2  = rbm.cascade_osd.avg_f2_ops;
        i = min([numel(lcc), numel(cllo), numel(cosd_f2m)]);
        ratio = 100 * (cllo(i) - lcc(i)) / max(lcc(i), 1);
        fprintf('  @ %.2f dB: 级联 LLOSD vs 纯 RS-LCC-BR 时延增幅 = %+.1f%%  (KPI <=+10%%: %s)\n', ...
            ebn0(i), ratio, ternary(ratio <= 10, '达标 OK', '未达标'));
        % Lagrange vs Gaussian-elim: compare TOTAL elementary ops (f2m + f2).
        % OSD's Gaussian elimination cost lives in f2, which LLOSD avoids.
        cllo_tot = cllo(i) + cllo_f2(i);
        cosd_tot = cosd_f2m(i) + cosd_f2(i);
        if cllo_tot < cosd_tot
            fprintf('  Lagrange 相对高斯消元省算力: LLOSD 总运算=%.0f < OSD 总运算=%.0f  (低时延成立)\n', ...
                cllo_tot, cosd_tot);
        else
            fprintf('  LLOSD 总运算=%.0f  OSD 总运算=%.0f\n', cllo_tot, cosd_tot);
        end
    end
    fprintf('====================================================\n');
end

function y = valAt(rbm, key, i, field)
    if nargin < 4, field = 'avg_f2m_ops'; end
    r = rbm.(key);
    v = r.(field);
    if i <= numel(v), y = v(i); else, y = NaN; end
end

function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end
