function paths = ensure_asset_dirs()
%ENSURE_ASSET_DIRS Make local Work 02 output directories if absent.
paths = cascade.work02_paths();
folders = {paths.assets, paths.data, paths.figures, paths.logs};
for ii = 1:numel(folders)
    if ~isfolder(folders{ii})
        mkdir(folders{ii});
    end
end
end
