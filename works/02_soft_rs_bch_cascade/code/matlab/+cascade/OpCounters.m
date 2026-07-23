classdef OpCounters < handle
    %OPCOUNTERS Lightweight accounting object shared by the decoders.
    properties
        f2 = 0
        f2m = 0
        fp = 0
        latency_us = 0
        n_tvs = 0
    end
    methods
        function add(obj, other)
            obj.f2 = obj.f2 + other.f2;
            obj.f2m = obj.f2m + other.f2m;
            obj.fp = obj.fp + other.fp;
            obj.latency_us = obj.latency_us + other.latency_us;
            obj.n_tvs = obj.n_tvs + other.n_tvs;
        end
        function out = to_struct(obj)
            out = struct('f2', obj.f2, 'f2m', obj.f2m, 'fp', obj.fp, ...
                'latency_us', obj.latency_us, 'n_tvs', obj.n_tvs);
        end
    end
end
