function outputs = generate_figure_05(dataRoot, outputDir)
%GENERATE_FIGURE_05 Reproduce Figure 5 from the archived summary table.

arguments
    dataRoot (1,1) string
    outputDir (1,1) string = string(fullfile(fileparts(mfilename('fullpath')), 'output'))
end

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

summaryFile = fullfile(dataRoot, 'syn', 'outputs', ...
    'synthetic_wave_range_autocorr_summary.csv');
if ~isfile(summaryFile)
    error('Figure 5 summary file not found: %s', summaryFile);
end
if ~isfolder(outputDir)
    mkdir(outputDir);
end

T = readtable(summaryFile, 'TextType', 'string');
requiredVariables = { ...
    'analysisStatus', 'suppressWaveAmpOutputs', 'lambda_m', ...
    'lambdaEstB_rangeWseAutocorr_m', 'A_m', 'Aest_m'};
missingVariables = setdiff(requiredVariables, T.Properties.VariableNames);
if ~isempty(missingVariables)
    error('Summary table is missing required variables: %s', ...
        strjoin(missingVariables, ', '));
end

okMask = string(T.analysisStatus) == "ok";
okMask = okMask & ~logical(T.suppressWaveAmpOutputs);

lambdaTrue = T.lambda_m;
lambdaEstimated = T.lambdaEstB_rangeWseAutocorr_m;
amplitudeTrue = T.A_m;
amplitudeEstimated = T.Aest_m;

validWavelength = okMask & isfinite(lambdaTrue) & ...
    isfinite(lambdaEstimated) & isfinite(amplitudeTrue);
validAmplitude = okMask & isfinite(amplitudeTrue) & ...
    isfinite(amplitudeEstimated) & isfinite(lambdaTrue);

if ~any(validWavelength) || ~any(validAmplitude)
    error('No finite, unsuppressed rows are available for Figure 5.');
end

figureWidthIn = 5.5;
figureHeightIn = 2.75;
fontName = 'Arial';
axesFontSize = 9;
labelFontSize = 10;
colorbarFontSize = 9;
colorbarLabelFontSize = 10;
panelLabelFontSize = 10;
axesLineWidth = 0.75;
markerArea = 36;

fig = figure('Color', 'w', 'Units', 'inches', ...
    'Position', [0.5 0.5 figureWidthIn figureHeightIn], ...
    'PaperUnits', 'inches', ...
    'PaperPosition', [0 0 figureWidthIn figureHeightIn], ...
    'PaperSize', [figureWidthIn figureHeightIn], ...
    'InvertHardcopy', 'off');
cleanupFigure = onCleanup(@() close_valid_figure(fig));

layout = tiledlayout(fig, 1, 2, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

% Panel (a): imposed versus estimated wavelength.
ax1 = nexttile(layout, 1);
hold(ax1, 'on');
grid(ax1, 'on');
box(ax1, 'on');
scatter(ax1, lambdaTrue(validWavelength), lambdaEstimated(validWavelength), ...
    markerArea, amplitudeTrue(validWavelength), 'o', 'filled', ...
    'MarkerFaceColor', 'flat', 'MarkerEdgeColor', [0.2 0.2 0.2], ...
    'LineWidth', 0.6);
colormap(ax1, coolwarm2(256));
amplitudeLimits = finite_limits(amplitudeTrue(validWavelength));
clim(ax1, amplitudeLimits);
cb1 = colorbar(ax1);
cb1.Label.String = '$A$ (m)';
cb1.Label.Interpreter = 'latex';
cb1.FontName = fontName;
cb1.FontSize = colorbarFontSize;
cb1.Label.FontName = fontName;
cb1.Label.FontSize = colorbarLabelFontSize;
cb1.Ticks = linspace(amplitudeLimits(1), amplitudeLimits(2), 5);
cb1.TickLabels = compose('%.2f', cb1.Ticks);

wavelengthLimits = equal_axis_limits([ ...
    lambdaTrue(validWavelength); lambdaEstimated(validWavelength)]);
plot(ax1, wavelengthLimits, wavelengthLimits, 'k--', 'LineWidth', 0.9);
xlim(ax1, wavelengthLimits);
ylim(ax1, wavelengthLimits);
wavelengthTicks = ceil(wavelengthLimits(1)):floor(wavelengthLimits(2));
xticks(ax1, wavelengthTicks);
yticks(ax1, wavelengthTicks);
xlabel(ax1, '$\lambda$ (m)', 'Interpreter', 'latex');
ylabel(ax1, '$\lambda_{\mathrm{est}}$ (m)', 'Interpreter', 'latex');
ax1.XLabel.FontName = fontName;
ax1.XLabel.FontSize = labelFontSize;
ax1.YLabel.FontName = fontName;
ax1.YLabel.FontSize = labelFontSize;
axis(ax1, 'square');
add_panel_label(ax1, '(a)');

% Panel (b): imposed versus estimated amplitude.
ax2 = nexttile(layout, 2);
hold(ax2, 'on');
grid(ax2, 'on');
box(ax2, 'on');
scatter(ax2, amplitudeTrue(validAmplitude), amplitudeEstimated(validAmplitude), ...
    markerArea, lambdaTrue(validAmplitude), 'o', 'filled', ...
    'MarkerFaceColor', 'flat', 'MarkerEdgeColor', [0.2 0.2 0.2], ...
    'LineWidth', 0.6);
colormap(ax2, coolwarm2(256));
wavelengthColorLimits = finite_limits(lambdaTrue(validAmplitude));
clim(ax2, wavelengthColorLimits);
cb2 = colorbar(ax2);
cb2.Label.String = '$\lambda$ (m)';
cb2.Label.Interpreter = 'latex';
cb2.FontName = fontName;
cb2.FontSize = colorbarFontSize;
cb2.Label.FontName = fontName;
cb2.Label.FontSize = colorbarLabelFontSize;

amplitudeAxisLimits = equal_axis_limits([ ...
    amplitudeTrue(validAmplitude); amplitudeEstimated(validAmplitude)]);
plot(ax2, amplitudeAxisLimits, amplitudeAxisLimits, 'k--', 'LineWidth', 0.9);
xlim(ax2, amplitudeAxisLimits);
ylim(ax2, amplitudeAxisLimits);
tickStep = 0.02;
tickIndices = ceil(amplitudeAxisLimits(1) / tickStep): ...
    floor(amplitudeAxisLimits(2) / tickStep);
amplitudeTicks = tickIndices * tickStep;
xticks(ax2, amplitudeTicks);
yticks(ax2, amplitudeTicks);
xtickformat(ax2, '%.2f');
ytickformat(ax2, '%.2f');
xlabel(ax2, '$A$ (m)', 'Interpreter', 'latex');
ylabel(ax2, '$A_{\mathrm{est}}$ (m)', 'Interpreter', 'latex');
ax2.XLabel.FontName = fontName;
ax2.XLabel.FontSize = labelFontSize;
ax2.YLabel.FontName = fontName;
ax2.YLabel.FontSize = labelFontSize;
axis(ax2, 'square');
add_panel_label(ax2, '(b)');

set([ax1 ax2], 'FontName', fontName, 'FontSize', axesFontSize, ...
    'LineWidth', axesLineWidth, 'TickDir', 'out', 'Layer', 'top', ...
    'GridAlpha', 0.12);
panelLabels = findall(fig, 'Tag', 'panelLabel');
set(panelLabels, 'FontName', fontName, 'FontSize', panelLabelFontSize);

pngFile = fullfile(outputDir, 'Figure5.png');
pdfFile = fullfile(outputDir, 'Figure5.pdf');
drawnow;
exportgraphics(fig, pngFile, 'Resolution', 600);
exportgraphics(fig, pdfFile, 'ContentType', 'vector');

outputs = struct( ...
    'pngFile', pngFile, ...
    'pdfFile', pdfFile, ...
    'summaryFile', summaryFile, ...
    'plottedWavelengthCases', nnz(validWavelength), ...
    'plottedAmplitudeCases', nnz(validAmplitude));
end

function limits = finite_limits(values)
values = values(isfinite(values));
minimum = min(values);
maximum = max(values);
if minimum == maximum
    padding = max(abs(minimum) * 0.05, 1e-6);
    limits = [minimum - padding, maximum + padding];
else
    limits = [minimum, maximum];
end
end

function limits = equal_axis_limits(values)
values = values(isfinite(values));
minimum = min(values);
maximum = max(values);
if minimum == maximum
    padding = max(abs(minimum) * 0.075, 1e-3);
else
    padding = 0.075 * (maximum - minimum);
end
limits = [minimum - padding, maximum + padding];
end

function add_panel_label(ax, label)
text(ax, 0.02, 0.97, label, 'Units', 'normalized', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'FontWeight', 'bold', 'Tag', 'panelLabel');
end

function close_valid_figure(fig)
if isgraphics(fig)
    close(fig);
end
end
