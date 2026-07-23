function write_json(filename, value)
%WRITE_JSON Save a MATLAB value as UTF-8 JSON, creating parent directory.
parent = fileparts(filename);
if ~isempty(parent) && ~isfolder(parent)
    mkdir(parent);
end
fid = fopen(filename, 'w', 'n', 'UTF-8');
if fid < 0
    error('cascade:IO:Open', 'Could not open %s for writing.', filename);
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s', jsonencode(value));
end
