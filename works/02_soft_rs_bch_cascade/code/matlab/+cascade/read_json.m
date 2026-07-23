function value = read_json(filename)
%READ_JSON Read JSON without relying on current directory.
if ~isfile(filename)
    error('cascade:IO:MissingFile', 'Required data file does not exist: %s', filename);
end
value = jsondecode(fileread(filename));
end
