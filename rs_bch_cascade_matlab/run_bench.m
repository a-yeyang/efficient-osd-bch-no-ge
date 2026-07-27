function res = run_bench(method_name, encoder, decoder, rate, ebn0_list, ...
                         n_info_symbols, m_bits, opts)
%RUN_BENCH  Monte Carlo BER/FER/latency SNR sweep (port of
%   cascade_src/simulate.py run_bench).
%
%   encoder : @(msg) -> coded bits
%   decoder : @(llr, counters) -> [msg_hat, info]   (info unused here)
%   rate    : effective code rate (for PAM4 sigma)
%   opts    : struct with fields seed, min_frame_errors, max_frames, verbose
%
%   Returns struct res with per-SNR arrays: ebn0_db, ser, ber, fer, n_frames,
%   n_frame_errors, avg_f2m_ops, avg_f2_ops, avg_latency_us.

    if nargin < 8, opts = struct(); end
    if ~isfield(opts, 'seed'),             opts.seed = 0;             end
    if ~isfield(opts, 'min_frame_errors'), opts.min_frame_errors = 15; end
    if ~isfield(opts, 'max_frames'),       opts.max_frames = 200;     end
    if ~isfield(opts, 'verbose'),          opts.verbose = true;       end

    res = struct('ebn0_db', [], 'ser', [], 'ber', [], 'fer', [], ...
        'n_frames', [], 'n_frame_errors', [], 'avg_f2m_ops', [], ...
        'avg_f2_ops', [], 'avg_latency_us', []);

    for ei = 1:numel(ebn0_list)
        ebn0 = ebn0_list(ei);
        rng = RandStream('mt19937ar', 'Seed', opts.seed + round(ebn0 * 100));
        n_frames = 0; n_frame_errors = 0; n_bit_errors = 0; n_symbol_errors = 0;
        sum_f2m = 0; sum_f2 = 0; sum_lat = 0;
        t_start = tic;

        while n_frames < opts.max_frames
            msg = double(randi(rng, [0, bitshift(1, m_bits) - 1], 1, n_info_symbols));
            coded = encoder(msg);
            llr = pam4.runChannel(coded, ebn0, rate, rng);

            counters = OpCounters();
            t_dec = tic;
            msg_hat = decoder(llr, counters);
            dec_us = toc(t_dec) * 1e6;

            n_frames = n_frames + 1;
            n_sym_err = sum(msg_hat(:).' ~= msg(:).');
            if n_sym_err > 0, n_frame_errors = n_frame_errors + 1; end
            n_symbol_errors = n_symbol_errors + n_sym_err;

            msg_bits = symbolsToBits(msg, m_bits);
            msg_hat_bits = symbolsToBits(msg_hat, m_bits);
            n_bit_errors = n_bit_errors + sum(msg_bits ~= msg_hat_bits);

            sum_f2m = sum_f2m + counters.f2m;
            sum_f2  = sum_f2  + counters.f2;
            sum_lat = sum_lat + dec_us;

            if n_frame_errors >= opts.min_frame_errors && n_frames >= 100
                break;
            end
        end

        elapsed = toc(t_start);
        ser = n_symbol_errors / max(1, n_frames * n_info_symbols);
        ber = n_bit_errors / max(1, n_frames * n_info_symbols * m_bits);
        fer = n_frame_errors / max(1, n_frames);
        avg_f2m = sum_f2m / max(1, n_frames);
        avg_f2  = sum_f2  / max(1, n_frames);
        avg_lat = sum_lat / max(1, n_frames);

        if opts.verbose
            fprintf('  %-26s @ %.2f dB: FER=%.3e BER=%.3e f2m=%.0f f2=%.0f lat=%.0fus %d frames %.1fs\n', ...
                method_name, ebn0, fer, ber, avg_f2m, avg_f2, avg_lat, n_frames, elapsed);
        end

        res.ebn0_db(end+1) = ebn0;
        res.ser(end+1) = ser;
        res.ber(end+1) = ber;
        res.fer(end+1) = fer;
        res.n_frames(end+1) = n_frames;
        res.n_frame_errors(end+1) = n_frame_errors;
        res.avg_f2m_ops(end+1) = avg_f2m;
        res.avg_f2_ops(end+1) = avg_f2;
        res.avg_latency_us(end+1) = avg_lat;

        if fer < 1e-6 && n_frame_errors < 3
            break;
        end
    end
end


function bits = symbolsToBits(syms, m)
    syms = double(syms(:).');
    n = numel(syms);
    bits = zeros(1, n*m);
    for i = 1:n
        for b = 0:m-1
            bits((i-1)*m + b + 1) = bitand(bitshift(syms(i), -b), 1);
        end
    end
end
