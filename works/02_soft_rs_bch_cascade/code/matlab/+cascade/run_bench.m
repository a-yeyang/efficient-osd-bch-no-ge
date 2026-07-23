function result = run_bench(method_name, codec_pack, ebn0_list, n_info_symbols, m_bits_per_symbol, varargin)
%RUN_BENCH Monte-Carlo BER/SER/FER runner for a configured codec pack.
%   Name/value options: seed, min_frame_errors, max_frames, verbose.
parser = inputParser();
parser.addParameter('seed', 0, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('min_frame_errors', 30, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('max_frames', 5000, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('verbose', true, @(x) islogical(x) || isnumeric(x));
parser.parse(varargin{:});
opts = parser.Results;

result = struct('ebn0_db', [], 'ser', [], 'ber', [], 'fer', [], ...
    'n_frames', [], 'n_frame_errors', [], 'avg_f2m_ops', [], 'avg_latency_us', []);
for snrIndex = 1:numel(ebn0_list)
    ebn0 = double(ebn0_list(snrIndex));
    rng(opts.seed + round(ebn0 * 100), 'twister');
    nFrames = 0;
    nFrameErrors = 0;
    nBitErrors = 0;
    nSymbolErrors = 0;
    sumF2m = 0;
    sumLatency = 0;
    totalTimer = tic;
    while nFrames < opts.max_frames
        message = randi([0, 2^m_bits_per_symbol - 1], 1, n_info_symbols);
        coded = codec_pack.encoder(message);
        llr = cascade.run_channel(coded, ebn0, codec_pack.effective_rate);
        decodeTimer = tic;
        counters = cascade.OpCounters();
        [messageHat, ~] = codec_pack.decoder(llr, counters);
        elapsedUs = toc(decodeTimer) * 1e6;
        nFrames = nFrames + 1;
        symbolErrors = sum(messageHat ~= message);
        if symbolErrors > 0
            nFrameErrors = nFrameErrors + 1;
        end
        nSymbolErrors = nSymbolErrors + symbolErrors;
        messageBits = cascade.symbols_to_bits(message, m_bits_per_symbol);
        decodedBits = cascade.symbols_to_bits(messageHat, m_bits_per_symbol);
        nBitErrors = nBitErrors + sum(messageBits ~= decodedBits);
        sumF2m = sumF2m + counters.f2m;
        sumLatency = sumLatency + elapsedUs;
        if nFrameErrors >= opts.min_frame_errors && nFrames >= 100
            break;
        end
    end
    elapsed = toc(totalTimer);
    result.ebn0_db(end + 1) = ebn0;
    result.ser(end + 1) = nSymbolErrors / max(1, nFrames * n_info_symbols);
    result.ber(end + 1) = nBitErrors / max(1, nFrames * n_info_symbols * m_bits_per_symbol);
    result.fer(end + 1) = nFrameErrors / max(1, nFrames);
    result.n_frames(end + 1) = nFrames;
    result.n_frame_errors(end + 1) = nFrameErrors;
    result.avg_f2m_ops(end + 1) = sumF2m / max(1, nFrames);
    result.avg_latency_us(end + 1) = sumLatency / max(1, nFrames);
    if opts.verbose
        fprintf('  %s @ %.2f dB: FER=%.3e, BER=%.3e, avg_f2m=%.0f, lat=%.0f us, %d frames, %.1fs\n', ...
            method_name, ebn0, result.fer(end), result.ber(end), result.avg_f2m_ops(end), ...
            result.avg_latency_us(end), nFrames, elapsed);
    end
    if result.fer(end) < 1e-6 && nFrameErrors < 3
        break;
    end
end
end
