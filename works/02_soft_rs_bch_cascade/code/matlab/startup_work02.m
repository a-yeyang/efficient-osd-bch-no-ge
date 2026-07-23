function paths = startup_work02()
%STARTUP_WORK02 Add the Work 02 MATLAB tree to the current MATLAB session.
root = fileparts(mfilename('fullpath'));
addpath(root);
addpath(fullfile(root, 'experiments'));
addpath(fullfile(root, 'tests'));
paths = cascade.work02_paths();
end
