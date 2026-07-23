function run_all_tests()
%RUN_ALL_TESTS Backward-compatible alias for the canonical MATLAB test runner.
%
% Kept for users who previously invoked ``run_all_tests``.  The authoritative
% runners live under each work's ``code/matlab/tests`` directory and are
% orchestrated by ``run_all_matlab_tests``.

    run_all_matlab_tests();
end
