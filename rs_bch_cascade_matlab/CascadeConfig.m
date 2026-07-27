classdef CascadeConfig
%CASCADECONFIG  Configuration for a cascaded RS+BCH scheme.
%
%   Port of cascade_src/cascade.py (CascadeConfig dataclass).

    properties
        m           % GF exponent (2^m - 1 = block length)
        k_rs        % RS message length (symbols)
        t_bch       % BCH error-correction capability
        llosd_tau = 2   % LLOSD decoding order for BCH inner
        lcc_eta   = 6   % LCC-BR eta for RS outer
    end

    methods
        function obj = CascadeConfig(m, k_rs, t_bch, llosd_tau, lcc_eta)
            obj.m = m;
            obj.k_rs = k_rs;
            obj.t_bch = t_bch;
            if nargin >= 4 && ~isempty(llosd_tau), obj.llosd_tau = llosd_tau; end
            if nargin >= 5 && ~isempty(lcc_eta), obj.lcc_eta = lcc_eta; end
        end

        function n = n_bch(obj)
            n = bitshift(1, obj.m) - 1;
        end

        function n = n_rs(obj)
            n = bitshift(1, obj.m) - 1;
        end

        function r = totalRate(obj)
            bch = BCHCode(obj.m, obj.t_bch);
            r = (obj.k_rs / obj.n_rs()) * (bch.k / bch.n);
        end

        function s = describe(obj)
            bch = BCHCode(obj.m, obj.t_bch);
            s = sprintf('RS(%d,%d) + BCH(%d,%d), total rate = %.4f', ...
                obj.n_rs(), obj.k_rs, bch.n, bch.k, obj.totalRate());
        end
    end
end
