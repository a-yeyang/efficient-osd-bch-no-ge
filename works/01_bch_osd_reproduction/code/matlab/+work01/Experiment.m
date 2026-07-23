classdef Experiment
    %EXPERIMENT Shared, deterministic helpers for the MATLAB plot/table entry points.

    methods (Static)
        function fers = run_fer(code, decoder, ebn0List, minErrors, maxFrames, seed, label)
            if nargin < 6, seed = 0; end
            if nargin < 7, label = ''; end
            rate = code.k/code.n;
            zero = zeros(1, code.n);
            x0 = work01.Core.bpsk_modulate(zero);
            fers = zeros(size(ebn0List));
            for ii = 1:numel(ebn0List)
                ebn0 = ebn0List(ii);
                rng(seed + round(ebn0*100), 'twister');
                sigma = work01.Core.sigma_from_ebn0(ebn0, rate);
                frames = 0; errors = 0; timer = tic;
                while frames < maxFrames
                    L = work01.Core.llr_from_y(work01.Core.awgn_channel(x0, sigma), sigma);
                    [cHat, ~] = decoder(code, L);
                    frames = frames + 1;
                    if ~isequal(cHat, zero), errors = errors + 1; end
                    if errors >= minErrors && frames >= 200, break; end
                end
                fers(ii) = errors/max(1, frames);
                fprintf('  %s @ %.1f dB: FER = %.2e (%d/%d) %.1fs\n', ...
                    label, ebn0, fers(ii), errors, frames, toc(timer));
            end
        end

        function values = avg_nbch(code, tau, ebn0List, nTrials, seed)
            if nargin < 5, seed = 0; end
            rate = code.k/code.n; zero = zeros(1,code.n);
            x0 = work01.Core.bpsk_modulate(zero); values = zeros(size(ebn0List));
            for ii = 1:numel(ebn0List)
                ebn0 = ebn0List(ii); rng(seed + round(ebn0*100), 'twister');
                sigma = work01.Core.sigma_from_ebn0(ebn0, rate); samples = zeros(1,nTrials);
                for trial = 1:nTrials
                    L = work01.Core.llr_from_y(work01.Core.awgn_channel(x0,sigma),sigma);
                    [~, s] = work01.Core.llosd_fast(code, L, tau, true, false);
                    samples(trial) = s.n_bch_candidates;
                end
                values(ii) = mean(samples);
                fprintf('  (%d,%d), tau=%d, %.1f dB: avg N_BCH = %.2f\n', ...
                    code.n, code.k, tau, ebn0, values(ii));
            end
        end

        function value = avg_stat(code, decoder, ebn0List, nTrials, statName, seed, label)
            if nargin < 6, seed = 0; end
            if nargin < 7, label = ''; end
            rate = code.k/code.n; zero = zeros(1,code.n);
            x0 = work01.Core.bpsk_modulate(zero); value = zeros(size(ebn0List));
            for ii = 1:numel(ebn0List)
                ebn0 = ebn0List(ii); rng(seed + round(ebn0*100), 'twister');
                sigma = work01.Core.sigma_from_ebn0(ebn0, rate); samples = zeros(1,nTrials);
                for trial = 1:nTrials
                    L = work01.Core.llr_from_y(work01.Core.awgn_channel(x0,sigma),sigma);
                    [~, s] = decoder(code,L);
                    samples(trial) = work01.Experiment.get_stat(s, statName);
                end
                value(ii) = mean(samples);
                fprintf('  %s @ %.1f dB: avg %s = %.3f\n', label, ebn0, statName, value(ii));
            end
        end

        function [fer, avgOps] = measure_fer(code, decoder, ebn0, maxFrames, minErrors, seed)
            if nargin < 4, maxFrames = 800; end
            if nargin < 5, minErrors = 20; end
            if nargin < 6, seed = 1234; end
            rate = code.k/code.n; zero = zeros(1,code.n); x0 = work01.Core.bpsk_modulate(zero);
            rng(seed + round(ebn0*100), 'twister'); sigma = work01.Core.sigma_from_ebn0(ebn0,rate);
            errors = 0; frames = 0; ops = 0;
            while frames < maxFrames
                L = work01.Core.llr_from_y(work01.Core.awgn_channel(x0,sigma),sigma);
                [cHat,s] = decoder(code,L); frames = frames + 1;
                if ~isequal(cHat,zero), errors = errors + 1; end
                ops = ops + s.counters.f2 + s.counters.f2m;
                if errors >= minErrors && frames >= 200, break; end
            end
            fer = errors/max(1,frames); avgOps = ops/max(1,frames);
        end

        function [snr, ops] = find_snr_for_fer(code, decoder, targetFer, snrRange, maxFrames, minErrors)
            if nargin < 4, snrRange = [2,7]; end
            if nargin < 5, maxFrames = 800; end
            if nargin < 6, minErrors = 20; end
            snrs = snrRange(1):0.5:snrRange(2);
            for snr = snrs
                [fer,ops] = work01.Experiment.measure_fer(code,decoder,snr,maxFrames,minErrors,1234);
                if fer <= targetFer, return; end
            end
            snr = snrRange(2); ops = NaN;
        end

        function metrics = instrument(code, decoder, ebn0, nTrials, seed)
            if nargin < 5, seed = 0; end
            rate = code.k/code.n; x0 = ones(1,code.n);
            rng(seed + round(ebn0*100), 'twister'); sigma = work01.Core.sigma_from_ebn0(ebn0,rate);
            totals = [0,0,0,0];
            for ii = 1:nTrials
                L = work01.Core.llr_from_y(work01.Core.awgn_channel(x0,sigma),sigma);
                [~,s] = decoder(code,L); c = s.counters;
                totals = totals + [c.f2,c.f2m,c.fp,c.latency_us];
            end
            totals = totals/nTrials;
            metrics = struct('F2',totals(1),'F2m',totals(2),'FP',totals(3), ...
                'latency_us',totals(4));
        end

        function value = get_stat(s, name)
            if isfield(s,name), value = s.(name); return; end
            if strcmp(name,'n_teps') && isfield(s,'n_teps_llosd')
                value = s.n_teps_llosd;
            elseif strcmp(name,'n_bch_candidates') && isfield(s,'n_bch_llosd')
                value = s.n_bch_llosd;
            else
                value = 0;
            end
        end

        function data = map_from_pairs(names, values)
            data = containers.Map(names, values);
        end

        function key = json_field(name)
            % JSON-decode applies this same valid-name transform in MATLAB.
            key = matlab.lang.makeValidName(name);
        end
    end
end
