function outputs = generate_figure_07(dataRoot, outputDir)
%GENERATE_FIGURE_07 Reproduce Figure 7 from the associated dataset.
%
% Inputs:
%   dataRoot  - root folder containing the input data
%   outputDir - folder where generated figures are saved
%
% Example:
%   generate_figure_07("C:\path\to\figure_07_data", ...
%                     "C:\path\to\figure_07_output")

arguments
    dataRoot (1,1) string
    outputDir (1,1) string = string(fullfile( ...
        fileparts(mfilename('fullpath')), 'output'))
end

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

if ~isfolder(dataRoot)
    error('Input data folder does not exist: %s', dataRoot);
end

if ~isfolder(outputDir)
    mkdir(outputDir);
end

opts = struct( ...
    'dataRoot', dataRoot, ...
    'saveFigure', true, ...
    'outFigureFile', char(fullfile(outputDir, 'Figure7')));

outputs = plot_figure_7(opts);

end
