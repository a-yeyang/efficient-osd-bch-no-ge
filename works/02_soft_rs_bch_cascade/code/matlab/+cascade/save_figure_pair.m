function save_figure_pair(fig, base_filename)
%SAVE_FIGURE_PAIR Save a PNG and a vector PDF under a caller-supplied path.
parent = fileparts(base_filename);
if ~isfolder(parent)
    mkdir(parent);
end
exportgraphics(fig, [base_filename, '.png'], 'Resolution', 140);
exportgraphics(fig, [base_filename, '.pdf'], 'ContentType', 'vector');
end
