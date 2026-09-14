function outputs = generate_figure_04(dataRoot, outputDir)
%GENERATE_FIGURE_04 Reproduce Figure 4 from the archived synthetic case.

arguments
    dataRoot (1,1) string
    outputDir (1,1) string = string(fullfile(fileparts(mfilename('fullpath')), 'output'))
end

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

checkpointFile = fullfile(dataRoot, 'syn', 'outputs', ...
    'syn_solver_inputs_0pt88_case_amp0.02_w1.52_ckpt.mat');
if ~isfile(checkpointFile)
    error('Figure 4 checkpoint file not found: %s', checkpointFile);
end
if ~isfolder(outputDir)
    mkdir(outputDir);
end

outputBase = fullfile(outputDir, 'Figure4');
opts = struct('outFigureFile', char(outputBase), 'saveFigure', true);

fig = KLT_plot_refined_wave_solution_three_panel_final( ...
    checkpointFile, [], [], opts);
close(fig);

outputs = struct( ...
    'pngFile', outputBase + ".png", ...
    'pdfFile', outputBase + ".pdf", ...
    'checkpointFile', checkpointFile);
end
