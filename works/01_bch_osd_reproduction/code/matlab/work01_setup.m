function paths = work01_setup()
%WORK01_SETUP Add the MATLAB reproduction to path and create output folders.

paths = work01_paths();
addpath(genpath(paths.matlab_root));
for folder = {paths.assets_root, paths.data_dir, paths.figures_dir, paths.tables_dir}
    if ~isfolder(folder{1})
        mkdir(folder{1});
    end
end
end
