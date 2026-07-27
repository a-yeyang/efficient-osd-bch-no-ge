classdef pam4
%PAM4  PAM4 Gray modulation + AWGN + per-bit LLR (port of cascade_src/pam4.py).
%
%   Constellation (Gray, 2 bits/symbol, b1=MSB b0=LSB):
%       level | b1 b0
%       -3    | 0 0
%       -1    | 0 1
%       +1    | 1 1
%       +3    | 1 0
%   Average symbol energy E_s = (9+1+1+9)/4 = 5.
%
%   All methods are static; call as pam4.bitsToPam4(bits) etc.

    properties (Constant)
        E_S_AVG = 5.0;
        ALL_LEVELS = [-3.0, -1.0, +1.0, +3.0];
    end

    methods (Static)
        function out = bitsToPam4(bits)
            % Map binary vector (length 2N) to N PAM4 symbols.
            bits = double(bits(:).');
            assert(mod(numel(bits), 2) == 0);
            n_syms = numel(bits) / 2;
            out = zeros(1, n_syms);
            for i = 1:n_syms
                b1 = bits(2*i - 1);
                b0 = bits(2*i);
                out(i) = pam4.bitPairToLevel(b1, b0);
            end
        end

        function lvl = bitPairToLevel(b1, b0)
            if b1 == 0 && b0 == 0,      lvl = -3.0;
            elseif b1 == 0 && b0 == 1,  lvl = -1.0;
            elseif b1 == 1 && b0 == 1,  lvl = +1.0;
            else,                       lvl = +3.0;   % (1,0)
            end
        end

        function out = pam4ToBitsHard(y)
            % Hard-decision: nearest level -> bit pair.
            y = y(:).';
            out = zeros(1, 2*numel(y));
            lv = pam4.ALL_LEVELS;
            for i = 1:numel(y)
                [~, mi] = min(abs(lv - y(i)));
                level = lv(mi);
                switch level
                    case -3.0, b1 = 0; b0 = 0;
                    case -1.0, b1 = 0; b0 = 1;
                    case +1.0, b1 = 1; b0 = 1;
                    otherwise, b1 = 1; b0 = 0;  % +3
                end
                out(2*i - 1) = b1;
                out(2*i)     = b0;
            end
        end

        function s = sigmaFromEbn0(ebn0_db, rate)
            % sigma^2 = E_s * rate / (4 * Eb/N0).
            ebn0_lin = 10.0 ^ (ebn0_db / 10.0);
            sigma2 = pam4.E_S_AVG * rate / (4.0 * ebn0_lin);
            s = sqrt(sigma2);
        end

        function y = awgnChannel(x, sigma, rng)
            % rng is a MATLAB RandStream for reproducibility.
            y = x + sigma * randn(rng, size(x));
        end

        function out = bitLlr(y, sigma)
            % Per-bit LLR: out(2i-1)=LLR(b1|y_i), out(2i)=LLR(b0|y_i).
            % LLR = log P(b=0|y)/P(b=1|y). Positive => bit 0 more likely.
            y = y(:).';
            inv_2s2 = 1.0 / (2.0 * sigma * sigma);
            nn = numel(y);
            out = zeros(1, 2*nn);
            for i = 1:nn
                d_m3 = -inv_2s2 * (y(i) - (-3.0))^2;
                d_m1 = -inv_2s2 * (y(i) - (-1.0))^2;
                d_p1 = -inv_2s2 * (y(i) - (+1.0))^2;
                d_p3 = -inv_2s2 * (y(i) - (+3.0))^2;
                llr_b1 = pam4.lse2(d_m3, d_m1) - pam4.lse2(d_p1, d_p3);
                llr_b0 = pam4.lse2(d_m3, d_p3) - pam4.lse2(d_m1, d_p1);
                out(2*i - 1) = llr_b1;
                out(2*i)     = llr_b0;
            end
        end

        function v = lse2(a, b)
            % log(exp(a)+exp(b)), numerically stable.
            if a > b, mx = a; else, mx = b; end
            v = mx + log1p(exp(-abs(a - b)));
        end

        function llr = runChannel(bits, ebn0_db, rate, rng)
            % Modulate PAM4 -> AWGN -> per-bit LLR. Pads one bit if odd.
            bits = double(bits(:).');
            sigma = pam4.sigmaFromEbn0(ebn0_db, rate);
            n_orig = numel(bits);
            if mod(n_orig, 2) ~= 0
                bits_p = [bits, 0];
            else
                bits_p = bits;
            end
            x = pam4.bitsToPam4(bits_p);
            y = pam4.awgnChannel(x, sigma, rng);
            llr = pam4.bitLlr(y, sigma);
            llr = llr(1:n_orig);
        end
    end
end
