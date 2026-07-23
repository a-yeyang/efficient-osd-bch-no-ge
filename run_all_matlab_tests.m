function run_all_matlab_tests()
%RUN_ALL_MATLAB_TESTS Run the deterministic MATLAB validation for all works.
%
% Each work owns its own MATLAB implementation and test runner.  Resetting the
% path between runners makes the dependency boundary explicit and prevents a
% same-named helper in a later work from silently shadowing an earlier one.

    repository_root = fileparts(mfilename('fullpath'));
    test_specs = {
        '01_bch_osd_reproduction', 'run_work01_tests';
        '02_soft_rs_bch_cascade',  'run_work02_tests';
        '03_hard_rs_bch_cascade',  'run_work03_tests';
    };

    fprintf('=== BCH/RS-BCH MATLAB verification ===\n');
    for i = 1:size(test_specs, 1)
        work_name = test_specs{i, 1};
        runner_name = test_specs{i, 2};
        matlab_root = fullfile(repository_root, 'works', work_name, 'code', 'matlab');
        restoredefaultpath();
        addpath(genpath(matlab_root));
        fprintf('\n--- %s ---\n', work_name);
        assert(exist(runner_name, 'file') == 2, ...
            'Missing MATLAB test runner: %s', runner_name);
        feval(runner_name);
    end
    fprintf('\n=== All MATLAB work tests passed ===\n');
end
