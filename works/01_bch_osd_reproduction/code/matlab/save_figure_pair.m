function save_figure_pair(fig, baseFilename)
%SAVE_FIGURE_PAIR Save a MATLAB figure as 140-dpi PNG and vector PDF.

parent = fileparts(baseFilename);
if ~isfolder(parent), mkdir(parent); end
exportgraphics(fig, [baseFilename '.png'], 'Resolution', 140);
exportgraphics(fig, [baseFilename '.pdf'], 'ContentType', 'vector');
end
