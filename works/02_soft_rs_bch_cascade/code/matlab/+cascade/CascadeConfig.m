classdef CascadeConfig
    %CASCADECONFIG Parameters shared by the RS+BCH cascade variants.
    properties
        m
        k_rs
        t_bch
        llosd_tau = 2
        lcc_eta = 6
    end
    properties (Dependent)
        n_bch
        n_rs
        total_rate
    end
    methods
        function obj = CascadeConfig(m, k_rs, t_bch, varargin)
            validateattributes(m, {'numeric'}, {'scalar', 'integer', '>=', 1, '<=', 8});
            obj.m = double(m);
            obj.k_rs = double(k_rs);
            obj.t_bch = double(t_bch);
            if mod(numel(varargin), 2) ~= 0
                error('cascade:Config:NameValue', 'Optional arguments must be name/value pairs.');
            end
            for ii = 1:2:numel(varargin)
                key = char(varargin{ii});
                value = varargin{ii + 1};
                switch lower(key)
                    case {'llosd_tau', 'tau'}
                        obj.llosd_tau = double(value);
                    case {'lcc_eta', 'eta'}
                        obj.lcc_eta = double(value);
                    otherwise
                        error('cascade:Config:UnknownOption', 'Unknown CascadeConfig option: %s', key);
                end
            end
            assert(obj.k_rs >= 1 && obj.k_rs <= obj.n_rs, 'cascade:Config:Krs', 'k_rs must be in 1..n_rs.');
        end

        function value = get.n_bch(obj)
            value = 2^obj.m - 1;
        end
        function value = get.n_rs(obj)
            value = 2^obj.m - 1;
        end
        function value = get.total_rate(obj)
            bch = cascade.BCHCode(obj.m, obj.t_bch);
            value = (obj.k_rs / obj.n_rs) * (bch.k / bch.n);
        end
        function text = describe(obj)
            bch = cascade.BCHCode(obj.m, obj.t_bch);
            text = sprintf('RS(%d,%d) + BCH(%d,%d), total rate = %.4f', ...
                obj.n_rs, obj.k_rs, bch.n, bch.k, obj.total_rate);
        end
    end
end
