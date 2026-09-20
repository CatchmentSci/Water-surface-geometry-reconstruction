function outputs = generate_figure_03(outputDir)
%GENERATE_FIGURE_KH_VS_F Reproduce the accepted kh-Fr figure.

arguments
    outputDir (1,1) string = string(fullfile( ...
        fileparts(mfilename('fullpath')), 'output'))
end

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

if ~isfolder(outputDir)
    mkdir(outputDir);
end

outputBase = fullfile(outputDir,'Figure3');

opts = struct( ...
    'saveFigure',true, ...
    'outFigureFile',char(outputBase));

fig = plot_kh_vs_F(opts);
close(fig);

outputs = struct( ...
    'pngFile',outputBase + ".png", ...
    'pdfFile',outputBase + ".pdf");
end
