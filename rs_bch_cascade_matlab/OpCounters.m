classdef OpCounters < handle
%OPCOUNTERS  Mutable operation counters (port of src/osd.py OpCounters).
%
%   Handle class so decoders can accumulate into it in place, exactly like the
%   Python dataclass passed by reference.
%
%     f2         binary XOR / AND ops
%     f2m        GF(2^m) multiplications / additions
%     fp         floating-point ops (adds and comparisons)
%     latency_us measured decode latency (microseconds, tic/toc)
%     n_tvs      number of Chase test-vectors processed (LCC-BR)

    properties
        f2         = 0
        f2m        = 0
        fp         = 0
        latency_us = 0
        n_tvs      = 0
    end

    methods
        function reset(obj)
            obj.f2 = 0; obj.f2m = 0; obj.fp = 0;
            obj.latency_us = 0; obj.n_tvs = 0;
        end
    end
end
