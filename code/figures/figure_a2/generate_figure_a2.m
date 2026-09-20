function outputs = generate_figure_a2(outputDir)
%GENERATE_FIGURE_SWW_VS_depth Reproduce the SWW-depth relationship figure.

arguments
    outputDir (1,1) string = string(fullfile( ...
        fileparts(mfilename('fullpath')), 'output'))
end

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

if ~isfolder(outputDir)
    mkdir(outputDir);
end

outputBase = fullfile(outputDir,'FigureA2');

opts = struct( ...
    'saveFigure',true, ...
    'outFigureFile',char(outputBase));

fig = plot_SWW_vs_depth(opts);
close(fig);

outputs = struct( ...
    'pngFile',outputBase + ".png", ...
    'pdfFile',outputBase + ".pdf");
end
