function write_json_file(filename, value)
%WRITE_JSON_FILE Write UTF-8 JSON, creating the containing folder as needed.

parent = fileparts(filename);
if ~isfolder(parent), mkdir(parent); end
fid = fopen(filename, 'w', 'n', 'UTF-8');
if fid < 0
    error('work01:IO:OpenFailed', 'Could not open %s for writing.', filename);
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fwrite(fid, jsonencode(value), 'char');
fwrite(fid, newline, 'char');
end
