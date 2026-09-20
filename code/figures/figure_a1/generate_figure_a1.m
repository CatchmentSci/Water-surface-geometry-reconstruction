function outputs = generate_figure_a1(outputDir)
%GENERATE_FIGURE_Us_VS_SWW Reproduce the Us-SWW relationship figure.

arguments
    outputDir (1,1) string = string(fullfile( ...
        fileparts(mfilename('fullpath')), 'output'))
end

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

if ~isfolder(outputDir)
    mkdir(outputDir);
end

outputBase = fullfile(outputDir,'FigureA1');

opts = struct( ...
    'saveFigure',true, ...
    'outFigureFile',char(outputBase));

fig = plot_Us_vs_SWW(opts);
close(fig);

outputs = struct( ...
    'pngFile',outputBase + ".png", ...
    'pdfFile',outputBase + ".pdf");
end
