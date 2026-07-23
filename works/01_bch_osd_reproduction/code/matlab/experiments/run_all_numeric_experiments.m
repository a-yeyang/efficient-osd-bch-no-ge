function run_all_numeric_experiments(options)
%RUN_ALL_NUMERIC_EXPERIMENTS Execute every numerical MATLAB experiment in order.
%
% Full defaults match the Python reproduction and are computationally heavy.
% To make a bounded exploratory run, pass e.g. struct('max_frames_scale',.01,
% 'n_trials_scale',.01, 'n_trials_63',2, 'n_trials_31',2).

if nargin < 1, options = struct(); end
fig02_nbch(options);
fig03_04_fer(options);
fig05_comparison(options);
fig07_09_longcode(options);
fig07_09_hsd_rerun(options);
fig08_rate_sweep(options);
fig10_11_ops(options);
tables(options);
latency_analysis(options);
end
