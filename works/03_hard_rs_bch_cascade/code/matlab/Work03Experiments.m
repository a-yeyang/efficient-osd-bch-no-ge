classdef Work03Experiments
%WORK03EXPERIMENTS MATLAB entry points corresponding to main.py/main_v2.py.
%
%   Work03Experiments.main(true)      % bounded profile Monte-Carlo + assets
%   Work03Experiments.main(false)     % Python-sized experiment schedule
%   Work03Experiments.main_v2(true)   % v2 latency-only experiment + assets
%
% The optional second argument writeOutputs is useful for tests.  All files,
% when requested, are written under works/03_hard_rs_bch_cascade/assets.

methods(Static)
    function out = main(profile, writeOutputs)
        %MAIN Port of code/python/experiments/main.py.
        if nargin < 1 || isempty(profile)
            profile = true;
        end
        if nargin < 2 || isempty(writeOutputs)
            writeOutputs = true;
        end
        profile = logical(profile);
        writeOutputs = logical(writeOutputs);

        configs = Work03Experiments.configs();
        names = fieldnames(configs);
        allResults = struct();
        allLatency = struct();
        for i = 1:numel(names)
            name = names{i};
            cfg = configs.(name);
            fprintf('\n\n########## %s: %s ##########\n', name, cfg.description);
            allResults.(name) = Work03Experiments.run_experiment(name, cfg, profile);
            allLatency.(name) = HardCascade.latency_summary_v1(cfg);
        end
        out = struct('results', allResults, 'latency', allLatency, ...
            'profile', profile, 'n255', allResults.n255, 'n127', allResults.n127);
        if writeOutputs
            Work03Experiments.write_main_outputs(out);
        end
    end

    function out = main_v2(profile, writeOutputs)
        %MAIN_V2 Port of code/python/experiments/main_v2.py.
        %
        % Python v2 recomputes only latency results and uses the v1 BER data
        % as background context.  The MATLAB port therefore does not launch
        % a second Monte-Carlo run here.
        if nargin < 1 || isempty(profile)
            profile = true;
        end
        if nargin < 2 || isempty(writeOutputs)
            writeOutputs = true;
        end
        configs = Work03Experiments.configs();
        names = fieldnames(configs);
        allLatency = struct();
        for i = 1:numel(names)
            allLatency.(names{i}) = HardCascade.latency_summary_v2(configs.(names{i}));
        end
        Work03Experiments.print_v2_kpi_table(allLatency);
        out = struct('latency', allLatency, 'profile', logical(profile), ...
            'n255_cycles', allLatency.n255.direct_v2, ...
            'n127_cycles', allLatency.n127.direct_v2);
        if logical(writeOutputs)
            Work03Experiments.write_v2_outputs(out);
        end
    end

    function out = bench(m, kRs, profile)
        %BENCH Backward-compatible direct-cascade profile helper.
        if nargin < 3 || isempty(profile)
            profile = true;
        end
        cfg = HardCascade.hard_cascade_config(m, kRs);
        codec = HardCascade.cascade_create(m, kRs, 'direct');
        [snr, minErrors, maxFrames] = Work03Experiments.bench_schedule(m, logical(profile), false);
        out = HardCascade.run_bench(sprintf('bench/direct/%d', cfg.n_rs), codec, snr, ...
            kRs, m, 0, minErrors, maxFrames, true);
        out.cycles = HardCascade.latency_cycles(codec, 'v2');
    end

    function configs = configs()
        %CONFIGS The two exact HardCascadeConfig instances in main.py.
        configs = struct( ...
            'n255', HardCascade.hard_cascade_config(8, 239), ...
            'n127', HardCascade.hard_cascade_config(7, 113));
    end

    function codecs = build_codecs(cfg)
        %BUILD_CODECS Python main.py's pure/conv/direct method group.
        codecs = struct( ...
            'pure_rs', HardCascade.pure_rs_create(cfg.m, cfg.k_rs), ...
            'cascade_conv', HardCascade.cascade_create(cfg.m, cfg.k_rs, 'conv'), ...
            'cascade_direct', HardCascade.cascade_create(cfg.m, cfg.k_rs, 'direct'));
    end

    function results = run_experiment(configName, cfg, profile)
        %RUN_EXPERIMENT Execute all three hard-decision methods for one config.
        codecs = Work03Experiments.build_codecs(cfg);
        methods = fieldnames(codecs);
        results = struct();
        for i = 1:numel(methods)
            method = methods{i};
            codec = codecs.(method);
            [snr, minErrors, maxFrames] = Work03Experiments.bench_schedule( ...
                cfg.m, logical(profile), strcmp(method, 'pure_rs'));
            fprintf('\n--- %s / %s ---\n', configName, method);
            results.(method) = HardCascade.run_bench( ...
                sprintf('%s/%s', configName, method), codec, snr, cfg.k_rs, cfg.m, ...
                0, minErrors, maxFrames, true);
        end
    end

    function [snr, minErrors, maxFrames] = bench_schedule(m, profile, pureRs)
        %BENCH_SCHEDULE Full schedule mirrors main.py; profile is test-safe.
        if profile
            if m == 8
                snr = 8.0;
            else
                snr = 7.0;
            end
            minErrors = 1;
            maxFrames = 2;
            return;
        end
        if m == 8
            snr = 6.0:0.5:12.0;
        else
            snr = 5.0:0.5:11.5;
        end
        minErrors = 20;
        if pureRs
            maxFrames = 800;
        else
            maxFrames = 400;
        end
    end

    function write_main_outputs(out)
        %WRITE_MAIN_OUTPUTS Save data and figures under this work's assets.
        paths = HardCascade.ensure_asset_dirs();
        suffix = Work03Experiments.output_suffix(out.profile);
        names = {'n255', 'n127'};
        for i = 1:numel(names)
            name = names{i};
            cfg = Work03Experiments.configs();
            perConfig = struct('results', out.results.(name), ...
                'latency', out.latency.(name), 'config', cfg.(name));
            HardCascade.write_json(fullfile(paths.data_root, ...
                sprintf('matlab_%s_%s_results.json', suffix, name)), perConfig);
        end
        HardCascade.write_json(fullfile(paths.data_root, ...
            sprintf('matlab_%s_all_results.json', suffix)), ...
            struct('results', out.results, 'latency', out.latency));

        ferFigure = Work03Experiments.make_fer_plot(out.results);
        HardCascade.save_figure_pair(ferFigure, fullfile(paths.figures_root, ...
            sprintf('matlab_%s_fer_all', suffix)));
        close(ferFigure);
        latencyFigure = Work03Experiments.make_latency_bar(out.latency);
        HardCascade.save_figure_pair(latencyFigure, fullfile(paths.figures_root, ...
            sprintf('matlab_%s_latency_bars', suffix)));
        close(latencyFigure);
        fprintf('Saved Work 03 MATLAB main artifacts under %s\n', paths.assets_root);
    end

    function write_v2_outputs(out)
        %WRITE_V2_OUTPUTS Save v2 latency data and both Python-equivalent plots.
        paths = HardCascade.ensure_asset_dirs();
        suffix = Work03Experiments.output_suffix(out.profile);
        HardCascade.write_json(fullfile(paths.data_root, ...
            sprintf('matlab_%s_v2_latency.json', suffix)), out.latency);
        bars = Work03Experiments.make_v2_latency_plot(out.latency);
        HardCascade.save_figure_pair(bars, fullfile(paths.figures_root, ...
            sprintf('matlab_%s_latency_bars_v2', suffix)));
        close(bars);
        evolution = Work03Experiments.make_evolution_plot(out.latency);
        HardCascade.save_figure_pair(evolution, fullfile(paths.figures_root, ...
            sprintf('matlab_%s_latency_evolution', suffix)));
        close(evolution);
        fprintf('Saved Work 03 MATLAB v2 artifacts under %s\n', paths.assets_root);
    end

    function suffix = output_suffix(profile)
        if profile
            suffix = 'profile';
        else
            suffix = 'full';
        end
    end

    function fig = make_fer_plot(allResults)
        %MAKE_FER_PLOT MATLAB equivalent of main.py's make_fer_plot.
        names = fieldnames(allResults);
        fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1200 450]);
        styles = {'pure_rs', 'd-', [0.85 0.1 0.1], 'Pure RS-BM'; ...
            'cascade_conv', 'o-', [0.1 0.25 0.9], 'Cascade (BCH-Conv + RS-BM)'; ...
            'cascade_direct', '^-', [0.1 0.55 0.2], 'Cascade (BCH-Direct + RS-BM)'};
        configs = Work03Experiments.configs();
        for panel = 1:numel(names)
            name = names{panel};
            subplot(1, numel(names), panel);
            hold on;
            for row = 1:size(styles, 1)
                method = styles{row, 1};
                result = allResults.(name).(method);
                semilogy(result.ebn0_db, max(result.fer, 1e-6), styles{row, 2}, ...
                    'Color', styles{row, 3}, 'DisplayName', styles{row, 4}, ...
                    'MarkerSize', 6, 'LineWidth', 1.1);
            end
            xlabel('E_b/N_0 (dB)');
            ylabel('FER');
            title(sprintf('%s: %s', name, configs.(name).description), 'Interpreter', 'none');
            grid on;
            legend('Location', 'best', 'FontSize', 8);
            hold off;
        end
    end

    function fig = make_latency_bar(allLatency)
        %MAKE_LATENCY_BAR MATLAB equivalent of main.py's latency bar chart.
        names = fieldnames(allLatency);
        fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1000 500]);
        axesHandle = axes(fig);
        hold(axesHandle, 'on');
        x = 1:numel(names);
        width = 0.16;
        metrics = {'pure_rs', 'Pure RS-BM (baseline)', [0.85 0.1 0.1]; ...
            'conv_no_share', 'Cascade Conv (no share)', [0.1 0.25 0.9]; ...
            'direct_no_share', 'Cascade Direct (no share)', [0.1 0.55 0.2]; ...
            'conv_lagrange', 'Cascade Conv (Lagrange shared)', [0.0 0.1 0.55]; ...
            'direct_lagrange', 'Cascade Direct (Lagrange shared)', [0.0 0.35 0.12]};
        for metric = 1:size(metrics, 1)
            values = zeros(1, numel(names));
            for i = 1:numel(names)
                values(i) = allLatency.(names{i}).(metrics{metric, 1});
            end
            bar(axesHandle, x + (metric - 3) * width, values, width, ...
                'DisplayName', metrics{metric, 2}, 'FaceColor', metrics{metric, 3}, ...
                'FaceAlpha', 0.85);
        end
        for i = 1:numel(names)
            target = allLatency.(names{i}).kpi_target;
            line(axesHandle, [x(i) - 2.5*width, x(i) + 2.5*width], [target target], ...
                'Color', 'k', 'LineStyle', '--', 'LineWidth', 1, 'HandleVisibility', 'off');
            text(axesHandle, x(i) + 0.05, target + 0.5, sprintf('+10%% KPI=%d', target), ...
                'FontSize', 8);
        end
        set(axesHandle, 'XTick', x, 'XTickLabel', names);
        ylabel(axesHandle, 'Clock cycles');
        title(axesHandle, 'Cascade latency: clock cycles vs Pure RS-BM baseline');
        grid(axesHandle, 'on');
        legend(axesHandle, 'Location', 'northwest', 'FontSize', 8);
        hold(axesHandle, 'off');
    end

    function fig = make_v2_latency_plot(allLatency)
        %MAKE_V2_LATENCY_PLOT Port of main_v2.py's bar plot.
        names = fieldnames(allLatency);
        fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1300 500]);
        for panel = 1:numel(names)
            name = names{panel};
            summary = allLatency.(name);
            subplot(1, numel(names), panel);
            labels = {'Pure RS-BM', 'Conv / no share', 'Direct / no share', ...
                'Direct / v1 share', 'Direct / v2 share'};
            values = [summary.pure_rs, summary.conv_none, summary.direct_none, ...
                summary.direct_v1, summary.direct_v2];
            colors = [0.85 0.1 0.1; 0.1 0.25 0.9; 0.56 0.93 0.56; ...
                0.24 0.70 0.44; 0.13 0.55 0.13];
            bars = bar(values, 'FaceColor', 'flat');
            bars.CData = colors;
            hold on;
            yline(summary.kpi_ceiling, 'k--', ...
                sprintf('+10%% KPI = %.1f cyc', summary.kpi_ceiling), 'LineWidth', 1);
            for i = 1:numel(values)
                if i == 1
                    label = 'baseline';
                else
                    label = sprintf('%+.1f%%', (values(i) - summary.pure_rs) / summary.pure_rs * 100);
                end
                text(i, values(i) + 0.35, sprintf('%d\n(%s)', values(i), label), ...
                    'HorizontalAlignment', 'center', 'FontSize', 8);
            end
            set(gca, 'XTick', 1:numel(labels), 'XTickLabel', labels, 'XTickLabelRotation', 18);
            ylabel('Clock cycles');
            title(sprintf('%s: BCH direct = %d cycles', name, summary.bch_direct_cyc));
            grid on;
            ylim([0, max(values) * 1.25]);
            hold off;
        end
        sgtitle('Cascade Latency v2: aggressive Lagrange sharing + smaller-GF Direct');
    end

    function fig = make_evolution_plot(allLatency)
        %MAKE_EVOLUTION_PLOT Port of main_v2.py's KPI-evolution plot.
        names = fieldnames(allLatency);
        fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1000 500]);
        axesHandle = axes(fig);
        stages = {'Conv / no share', 'Direct / no share', 'Direct / v1 share', 'Direct / v2 share'};
        hold(axesHandle, 'on');
        for i = 1:numel(names)
            summary = allLatency.(names{i});
            cycles = [summary.conv_none, summary.direct_none, summary.direct_v1, summary.direct_v2];
            ratios = (cycles - summary.pure_rs) / summary.pure_rs * 100;
            plot(axesHandle, 1:numel(stages), ratios, 'o-', 'LineWidth', 2, ...
                'MarkerSize', 8, 'DisplayName', sprintf('%s (RS baseline=%d cyc)', ...
                names{i}, summary.pure_rs));
            for point = 1:numel(ratios)
                text(axesHandle, point, ratios(point) + 1.4, sprintf('%+.1f%%', ratios(point)), ...
                    'HorizontalAlignment', 'center', 'FontSize', 8);
            end
        end
        yline(axesHandle, 10, 'r--', '10% KPI ceiling', 'LineWidth', 1.5);
        yline(axesHandle, 0, ':', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
        set(axesHandle, 'XTick', 1:numel(stages), 'XTickLabel', stages);
        ylabel(axesHandle, 'Latency increase vs Pure RS-BM (%)');
        title(axesHandle, 'Cascade latency evolution: Conv to Direct v2 share');
        grid(axesHandle, 'on');
        ylim(axesHandle, [-5 55]);
        legend(axesHandle, 'Location', 'northeast', 'FontSize', 9);
        hold(axesHandle, 'off');
    end

    function print_v2_kpi_table(allLatency)
        %PRINT_V2_KPI_TABLE Console equivalent of main_v2.py's KPI table.
        fprintf('%s\n', repmat('=', 1, 90));
        fprintf('KPI Table (Clock Cycles) - v1 vs v2 optimization\n');
        fprintf('%s\n', repmat('=', 1, 90));
        names = fieldnames(allLatency);
        rows = {'Cascade Conv (no share)', 'conv_none'; ...
            'Cascade Conv + v1 share', 'conv_v1'; ...
            'Cascade Conv + v2 share', 'conv_v2'; ...
            'Cascade Direct (no share)', 'direct_none'; ...
            'Cascade Direct + v1 share', 'direct_v1'; ...
            'Cascade Direct + v2 share', 'direct_v2'};
        for i = 1:numel(names)
            name = names{i};
            summary = allLatency.(name);
            fprintf('\n### %s (RS-BM baseline = %d cyc, 10%% KPI ceiling = %.1f cyc)\n', ...
                name, summary.pure_rs, summary.kpi_ceiling);
            fprintf('BCH-Direct cycles: %d  (BCH-Conv: %d)\n', ...
                summary.bch_direct_cyc, summary.bch_conv_cyc);
            fprintf('%-40s %8s %12s %8s\n', 'Config', 'Cycles', 'vs baseline', 'KPI');
            for row = 1:size(rows, 1)
                cycles = summary.(rows{row, 2});
                ratio = (cycles - summary.pure_rs) / summary.pure_rs * 100;
                if ratio <= 10
                    status = 'PASS';
                else
                    status = 'FAIL';
                end
                fprintf('%-40s %8d %+11.1f%% %8s\n', rows{row, 1}, cycles, ratio, status);
            end
        end
    end
end
end
