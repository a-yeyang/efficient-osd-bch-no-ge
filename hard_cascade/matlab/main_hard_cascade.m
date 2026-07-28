function main_hard_cascade(fast)
    % One-click hard-decision RS+BCH cascade simulation.
    %
    % Pipeline (per frame): construct/encode cascade codeword -> PAM4 modulate
    % -> AWGN channel -> hard-decision demodulate -> BCH inner decode
    % -> RS outer decode -> compare recovered message bits (BER).
    %
    % Covers both v3 configs (report_v3.tex):
    %   Config 1 (KP4): RS(544,514,t=15)/GF(2^10) + BCH(144,136,t=1)/GF(2^8)
    %   Config 2 (255): RS(255,239,t=8)/GF(2^8)   + BCH(255,239,t=2)/GF(2^8)
    %
    % Usage:
    %   >> main_hard_cascade         % fast smoke run (few points/frames)
    %   >> main_hard_cascade(false)  % full BER-SNR sweep (slower)
    if nargin < 1 || isempty(fast), fast = true; end
    here = fileparts(mfilename('fullpath'));
    addpath(here);

    rng(42);

    ok = smoke_test_v3();
    if ~ok
        error('main_hard_cascade: 正确性自检未通过，已中止仿真');
    end
    fprintf('\n');

    config_names = {'cfg1_kp4', 'cfg2_255'};
    if fast
        ebn0_list = 6:2:12;      % dB sweep
        n_frames = 20;           % Monte-Carlo frames per Eb/N0 point
    else
        ebn0_list = 4:1:12;
        n_frames = 200;
    end

    fprintf('=== 硬判决 RS+BCH 级联码：一键仿真 ===\n\n');

    results = struct();
    for ci = 1:numel(config_names)
        cname = config_names{ci};
        cfg = cascade_config(cname);
        codec = cascade_init(cfg, 'direct');

        fprintf('--- %s ---\n', cfg.name);
        fprintf('RS(%d,%d,t=%d)/GF(2^%d) + BCH(%d,%d,t=%d)/GF(2^%d)\n', ...
            cfg.n_rs, cfg.k_rs, cfg.t_rs, cfg.m_rs, ...
            cfg.n_bch, cfg.k_bch, cfg.t_bch, cfg.m_bch);
        fprintf('BCH blocks: %d, pad bits: %d, coded bits: %d, rate: %.4f\n', ...
            codec.n_bch_blocks, codec.n_pad_bits, codec.n_coded_bits, codec.effective_rate);

        ber = zeros(1, numel(ebn0_list));
        fer = zeros(1, numel(ebn0_list));
        for ei = 1:numel(ebn0_list)
            ebn0_db = ebn0_list(ei);
            n_bit_err = 0;
            n_bit_tot = 0;
            n_frame_err = 0;
            for f = 1:n_frames
                % 1. Construct + encode
                msg_symbols = randi([0, (2^cfg.m_rs)-1], 1, cfg.k_rs);
                coded_bits = cascade_encode(codec, msg_symbols);

                % 2-4. PAM4 modulate -> AWGN -> hard-decision demodulate
                hard_bits = pam4_channel_hard(coded_bits, ebn0_db, codec.effective_rate);

                % 5. BCH inner decode + RS outer decode
                [msg_hat, ok_rs, ~] = cascade_decode(codec, hard_bits); %#ok<ASGLU>

                tx_bits = symbols_to_bits_lsb(msg_symbols, cfg.m_rs);
                rx_bits = symbols_to_bits_lsb(msg_hat, cfg.m_rs);
                n_err = sum(tx_bits ~= rx_bits);
                n_bit_err = n_bit_err + n_err;
                n_bit_tot = n_bit_tot + numel(tx_bits);
                if ~ok_rs || n_err > 0
                    n_frame_err = n_frame_err + 1;
                end
            end
            ber(ei) = n_bit_err / n_bit_tot;
            fer(ei) = n_frame_err / n_frames;
            fprintf('  Eb/N0=%5.1f dB   BER=%.3e   FER=%.3e\n', ebn0_db, ber(ei), fer(ei));
        end
        results.(cname).ebn0 = ebn0_list;
        results.(cname).ber = ber;
        results.(cname).fer = fer;
        results.(cname).cfg = cfg;
        fprintf('\n');
    end

    out_dir = fullfile(here, 'results');
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
    save(fullfile(out_dir, 'main_hard_cascade_results.mat'), 'results');
    fprintf('结果已保存: %s\n', fullfile(out_dir, 'main_hard_cascade_results.mat'));

    try
        figure('Visible', 'off');
        semilogy(results.cfg1_kp4.ebn0, max(results.cfg1_kp4.ber, 1e-12), '-o', ...
                 results.cfg2_255.ebn0, max(results.cfg2_255.ber, 1e-12), '-s');
        grid on;
        xlabel('Eb/N0 (dB)'); ylabel('BER');
        legend('Config 1 (KP4)', 'Config 2 (255)', 'Location', 'southwest');
        title('硬判决 RS+BCH 级联码 BER-SNR（一键仿真）');
        saveas(gcf, fullfile(out_dir, 'ber_snr_hard_cascade.png'));
        fprintf('BER 曲线已保存: %s\n', fullfile(out_dir, 'ber_snr_hard_cascade.png'));
    catch err
        fprintf('绘图跳过（无图形环境）: %s\n', err.message);
    end

    fprintf('\n=== 仿真完成 ===\n');
end
