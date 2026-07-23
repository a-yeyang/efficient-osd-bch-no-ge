function paths = work02_paths()
%WORK02_PATHS Resolve Work 02 directories from this file, never from pwd.
thisFile = mfilename('fullpath');
packageDir = fileparts(thisFile);                  % .../code/matlab/+cascade
matlabDir = fileparts(packageDir);                 % .../code/matlab
codeDir = fileparts(matlabDir);                    % .../code
workRoot = fileparts(codeDir);                     % .../02_soft_rs_bch_cascade
paths = struct();
paths.matlab_root = matlabDir;
paths.work_root = workRoot;
paths.assets = fullfile(workRoot, 'assets');
paths.data = fullfile(paths.assets, 'data');
paths.figures = fullfile(paths.assets, 'figures');
paths.logs = fullfile(paths.assets, 'logs');
end
