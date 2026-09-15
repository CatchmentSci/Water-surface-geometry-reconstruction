function outputs = generate_figure_11(archiveRoot, outputFolder, options)
%GENERATE_FIGURE_11 Amplification of surface-velocity bias in depth inversion.
%   OUTPUTS = GENERATE_FIGURE_11(ARCHIVEROOT, OUTPUTFOLDER) reads the
%   bundled initial-planar-velocity inverse-depth summary and the Table C1
%   wavelengths from the extracted Zenodo archive, and writes a two-panel
%   publication PNG and PDF to OUTPUTFOLDER.
%
%   For a fractional underestimation epsilon of the tracked surface
%   velocity, the constant-profile inversion (Eq. 10a) returns a depth
%   h_est satisfying tanh(k h_est) = (1-epsilon)^2 tanh(k h), so that
%
%       h_est/h = atanh((1-epsilon)^2 * tanh(kh)) / kh .
%
%   Panel (a) plots this predicted ratio for several epsilon against
%   relative depth kh = 2*pi*h/lambda, together with the observed ratio for
%   the stable R1-R11 cases. h and h_est are the paired case medians used in
%   Figure 9(a); lambda is the accepted autocorrelation wavelength reported
%   in Table C1. Grey bars show the interquartile ranges of h (mapped to kh)
%   and of h_est (divided by the median h). Marker colour denotes discharge.
%
%   Panel (b) plots the same ratio as a function of epsilon for fixed
%   values of kh, with the reference bias marked.
%
%   Name-value options:
%     "EpsilonValues"   bias fractions drawn in panel (a)
%                       (default [0.08 0.12 0.16]; first is the solid line)
%     "KhValues"        relative depths drawn in panel (b)
%                       (default [1 3 8])
%     "ReferenceEpsilon" bias marked in panel (b). Empty (default) uses the
%                       median initial-planar-velocity bias from Figure 8a.

arguments
    archiveRoot (1, 1) string
    outputFolder (1, 1) string = fullfile(fileparts(mfilename('fullpath')), "output")
    options.EpsilonValues (1, :) double {mustBePositive, mustBeLessThan(options.EpsilonValues, 1)} = [0.08 0.12 0.16]
    options.KhValues (1, :) double {mustBePositive} = [1 3 8]
    options.ReferenceEpsilon (1, :) double = []
end

archiveRoot = string(archiveRoot);
outputFolder = string(outputFolder);
scriptFolder = string(fileparts(mfilename("fullpath")));
outputRoot = fullfile(archiveRoot, "videos", "outputs");
depthFile = fullfile(scriptFolder, "..", "figure_09", "data", ...
    "wse_autocorrelation_uniform_linear_power_depth_summary.csv");
velocityFile = fullfile(scriptFolder, "..", "figure_08", "data", ...
    "wse_autocorrelation_velocity_method_sensitivity_summary.mat");
tableC1File = fullfile(outputRoot, "table_c1_real_wave_hydraulic_results.csv");
requiredFiles = [depthFile, velocityFile, tableC1File];
for ii = 1:numel(requiredFiles)
    if ~isfile(requiredFiles(ii))
        error("Figure11:MissingInput", ...
            "Required archive file was not found:\n%s", requiredFiles(ii));
    end
end

depth = readtable(depthFile, "VariableNamingRule", "preserve");
tableC1 = readtable(tableC1File, "VariableNamingRule", "preserve");

assert_fields(depth, ["caseNumber", "caseLabel", "discharge_Q_m3ps", ...
    "includedInFigure", "uniform_deep_nPairs", "uniform_deep_xMedian", ...
    "uniform_deep_xQ25", "uniform_deep_xQ75", "uniform_deep_yMedian", ...
    "uniform_deep_yQ25", "uniform_deep_yQ75"], "depth summary");
assert_fields(tableC1, ["caseNumber", "caseLabel", "lambdaEst_m"], "Table C1");

if height(depth) ~= 13
    error("Figure11:UnexpectedCases", ...
        "Expected the complete 13-case R1-R13 depth summary.");
end
expectedLabels = "R" + string((1:13).');
if ~isequal(string(depth.caseLabel), expectedLabels)
    error("Figure11:UnexpectedCases", ...
        "Depth summary is not the expected ordered R1-R13 set.");
end
included = logical(depth.includedInFigure(:));
if nnz(included) ~= 11 || ~all(included(1:11)) || any(included(12:13))
    error("Figure11:UnexpectedExclusions", ...
        "Expected 11 plotted cases with R12 and R13 excluded.");
end
if height(tableC1) ~= 11 || ...
        ~isequal(string(tableC1.caseLabel), expectedLabels(1:11))
    error("Figure11:UnexpectedCases", ...
        "Table C1 is not the expected ordered R1-R11 set.");
end

% Align Table C1 (R1-R11) with the depth summary (R1-R13).
lambda = nan(13, 1);
lambda(1:11) = double(tableC1.lambdaEst_m);

q = double(depth.discharge_Q_m3ps);
nPairs = double(depth.uniform_deep_nPairs);
h = double(depth.uniform_deep_xMedian);
hQ25 = double(depth.uniform_deep_xQ25);
hQ75 = double(depth.uniform_deep_xQ75);
hEst = double(depth.uniform_deep_yMedian);
hEstQ25 = double(depth.uniform_deep_yQ25);
hEstQ75 = double(depth.uniform_deep_yQ75);

valid = included & nPairs > 0 & isfinite(lambda) & lambda > 0 & ...
    isfinite(h) & h > 0 & isfinite(hEst) & isfinite(q) & ...
    isfinite(hQ25) & isfinite(hQ75) & isfinite(hEstQ25) & isfinite(hEstQ75);

k = 2 .* pi ./ lambda;
kh = k .* h;
ratio = hEst ./ h;
% Interquartile ranges: kh from the surveyed-depth quartiles; the ratio
% from the inverse-depth quartiles divided by the median surveyed depth.
khErrLow = kh - k .* hQ25;
khErrHigh = k .* hQ75 - kh;
ratioErrLow = ratio - hEstQ25 ./ h;
ratioErrHigh = hEstQ75 ./ h - ratio;

epsilonValues = options.EpsilonValues(:).';
predictedAtCases = nan(13, numel(epsilonValues));
for ee = 1:numel(epsilonValues)
    predictedAtCases(:, ee) = predicted_ratio(kh, epsilonValues(ee));
end
% Velocity bias that would reproduce each observed ratio exactly.
impliedEpsilon = 1 - sqrt(tanh(ratio .* kh) ./ tanh(kh));

if isempty(options.ReferenceEpsilon)
    velocitySummary = load(velocityFile, ...
        "autocorrCaseSummaryTable", "autocorrPanelIncluded");
    velocityCases = velocitySummary.autocorrCaseSummaryTable;
    velocityIncluded = logical(velocitySummary.autocorrPanelIncluded(:));
    referenceValues = 1 - double(velocityCases.deep_xMedian) ./ ...
        double(velocityCases.deep_yMedian);
    referenceEpsilon = round(median(referenceValues(velocityIncluded & ...
        isfinite(referenceValues)), "omitnan"), 2);
else
    if ~isscalar(options.ReferenceEpsilon) || ...
            ~isfinite(options.ReferenceEpsilon) || ...
            options.ReferenceEpsilon < 0 || options.ReferenceEpsilon >= 1
        error("Figure11:InvalidReferenceEpsilon", ...
            "ReferenceEpsilon must be empty or a scalar in [0, 1).")
    end
    referenceEpsilon = options.ReferenceEpsilon;
end

figCfg = struct;
figCfg.figureWidth_in = 5.5;
figCfg.figureHeight_in = 3.15;
figCfg.figureResolution_dpi = 600;
figCfg.axesFontSize = 9;
figCfg.labelFontSize = 10;
figCfg.colorbarFontSize = 9;
figCfg.legendFontSize = 7.5;
figCfg.panelLabelFontSize = 10;
figCfg.annotationFontSize = 7.5;
figCfg.markerArea = 30;
figCfg.errorColour = [0.62 0.62 0.62];
figCfg.curveColour = [0.10 0.10 0.10];
figCfg.curveStyles = ["-", "--", ":", "-."];

finiteQ = q(valid);
if isempty(finiteQ)
    qLimits = [0 1];
else
    qLimits = [min(finiteQ), max(finiteQ)];
    if qLimits(1) == qLimits(2)
        qLimits = qLimits + [-0.5 0.5];
    end
end

if ~isfolder(outputFolder)
    mkdir(outputFolder);
end
outputBase = fullfile(outputFolder, "wse_autocorrelation_depth_bias_amplification");

fig = figure("Color", "w", "Units", "inches", ...
    "Position", [1 1 figCfg.figureWidth_in figCfg.figureHeight_in], ...
    "PaperPositionMode", "auto");
cleanup = onCleanup(@() close(fig));
layout = tiledlayout(fig, 1, 2, ...
    "TileSpacing", "compact", "Padding", "compact");
colormap(fig, local_coolwarm(256));

ax1 = nexttile(layout, 1);
legendHandle = plot_amplification_panel(ax1, kh(valid), ratio(valid), ...
    khErrLow(valid), khErrHigh(valid), ratioErrLow(valid), ...
    ratioErrHigh(valid), q(valid), qLimits, epsilonValues, figCfg);

ax2 = nexttile(layout, 2);
plot_sensitivity_panel(ax2, options.KhValues, referenceEpsilon, ...
    figCfg);

cb = colorbar(ax1, "southoutside");
cb.Layout.Tile = "south";
cb.Label.Interpreter = "latex";
cb.Label.String = "$Q$ ($\mathrm{m}^{3}\,\mathrm{s}^{-1}$)";
cb.Label.FontSize = figCfg.labelFontSize;
cb.FontSize = figCfg.colorbarFontSize;
cb.TickLabelInterpreter = "tex";
cb.TickDirection = "out";
position_curve_legend(ax1, legendHandle);

pngFile = outputBase + ".png";
pdfFile = outputBase + ".pdf";
csvFile = outputBase + ".csv";
exportgraphics(fig, pngFile, "Resolution", figCfg.figureResolution_dpi);
exportgraphics(fig, pdfFile, "ContentType", "vector");

% Deposit the plotted values beside the figure.
caseTable = table(double(depth.caseNumber), string(depth.caseLabel), ...
    q, lambda, h, hQ25, hQ75, hEst, hEstQ25, hEstQ75, kh, ratio, ...
    impliedEpsilon, double(valid), ...
    VariableNames=["caseNumber", "caseLabel", "discharge_Q_m3ps", ...
    "lambdaEst_m", "h_pairedMedian_m", "h_Q25_m", "h_Q75_m", ...
    "hEst_pairedMedian_m", "hEst_Q25_m", "hEst_Q75_m", "kh", ...
    "hEst_over_h", "impliedVelocityBias", "plotted"]);
for ee = 1:numel(epsilonValues)
    caseTable.(sprintf("predicted_eps%02d", round(100*epsilonValues(ee)))) = ...
        predictedAtCases(:, ee);
end
writetable(caseTable, csvFile);

outputs = struct;
outputs.depthFile = depthFile;
outputs.velocityFile = velocityFile;
outputs.tableC1File = tableC1File;
outputs.pngFile = pngFile;
outputs.pdfFile = pdfFile;
outputs.csvFile = csvFile;
outputs.epsilonValues = epsilonValues;
outputs.khValues = options.KhValues;
outputs.referenceEpsilon = referenceEpsilon;
outputs.caseTable = caseTable;
outputs.includedCases = string(depth.caseLabel(valid));
outputs.excludedCases = string(depth.caseLabel(~valid));
end

function r = predicted_ratio(kh, epsilon)
%PREDICTED_RATIO h_est/h for a fractional velocity underestimate epsilon.
%   Solves tanh(k h_est) = (1-epsilon)^2 tanh(k h) for h_est and divides by h.
r = atanh((1 - epsilon).^2 .* tanh(kh)) ./ kh;
end

function assert_fields(T, requiredFields, description)
missingFields = setdiff(requiredFields, string(T.Properties.VariableNames));
if ~isempty(missingFields)
    error("Figure11:InvalidInput", "%s is missing field(s): %s", ...
        description, strjoin(missingFields, ", "));
end
end

function legendHandle = plot_amplification_panel(ax, kh, ratio, ...
        khErrLow, khErrHigh, ratioErrLow, ratioErrHigh, q, qLimits, ...
        epsilonValues, cfg)
xLimits = padded_limits(kh - khErrLow, kh + khErrHigh);
xLimits(1) = 0;
yLimits = padded_limits(ratio - ratioErrLow, ratio + ratioErrHigh);
yLimits(1) = 0;
yLimits(2) = max(yLimits(2), 1.0);

hold(ax, "on");
grid(ax, "off");
box(ax, "on");
xlim(ax, xLimits);
ylim(ax, yLimits);

khCurve = linspace(max(0.05, 0.02 .* xLimits(2)), xLimits(2), 600);
legendHandles = gobjects(1, numel(epsilonValues));
legendLabels = strings(1, numel(epsilonValues));
for ee = 1:numel(epsilonValues)
    style = cfg.curveStyles(min(ee, numel(cfg.curveStyles)));
    lineWidth = 1.1;
    if ee == 1
        lineWidth = 1.5;
    end
    legendHandles(ee) = plot(ax, khCurve, ...
        predicted_ratio(khCurve, epsilonValues(ee)), style, ...
        "Color", cfg.curveColour, "LineWidth", lineWidth);
    legendLabels(ee) = sprintf("\\epsilon = %g%%", 100 .* epsilonValues(ee));
end

draw_xy_iqr(ax, kh, ratio, khErrLow, khErrHigh, ratioErrLow, ...
    ratioErrHigh, cfg.errorColour, yLimits);
scatter(ax, kh, ratio, cfg.markerArea, q, "filled", ...
    "MarkerEdgeColor", [0.15 0.15 0.15], "LineWidth", 0.5, ...
    "HandleVisibility", "off");
clim(ax, qLimits);

legendHandle = legend(ax, legendHandles, cellstr(legendLabels), ...
    "Location", "none", "FontSize", cfg.legendFontSize, "Box", "off", ...
    "Interpreter", "tex");
title(legendHandle, "Predicted", "FontSize", cfg.legendFontSize, ...
    "FontWeight", "normal");

xlabel(ax, "$kh=2\pi h/\lambda$", ...
    "Interpreter", "latex", "FontSize", cfg.labelFontSize);
ylabel(ax, "$h_{\mathrm{est}}/h$", ...
    "Interpreter", "latex", "FontSize", cfg.labelFontSize);
add_panel_label(ax, "(a)", cfg.panelLabelFontSize, [0.03 0.97], "top");
format_axes(ax, cfg);
end

function plot_sensitivity_panel(ax, khValues, referenceEpsilon, cfg)
epsilonCurve = linspace(0, 0.15, 400);
hold(ax, "on");
grid(ax, "off");
box(ax, "on");
xlim(ax, [0 15]);
ylim(ax, [0 1.05]);
if referenceEpsilon > 0
    plot(ax, 100 .* [referenceEpsilon referenceEpsilon], [0 1.05], ":", ...
        "Color", cfg.errorColour, "LineWidth", 0.8, "HandleVisibility", "off");
    text(ax, 100 .* referenceEpsilon - 0.4, 1.02, ...
        sprintf("%g%%", 100 .* referenceEpsilon), ...
        "FontSize", cfg.annotationFontSize, "Color", [0.35 0.35 0.35], ...
        "HorizontalAlignment", "right", "VerticalAlignment", "top", ...
        "Interpreter", "tex");
end
for kk = 1:numel(khValues)
    style = cfg.curveStyles(min(kk, numel(cfg.curveStyles)));
    r = predicted_ratio(khValues(kk), epsilonCurve);
    plot(ax, 100 .* epsilonCurve, r, style, "Color", cfg.curveColour, ...
        "LineWidth", 1.1, "HandleVisibility", "off");
    labelEpsilon = 0.78 .* epsilonCurve(end);
    labelRatio = predicted_ratio(khValues(kk), labelEpsilon) + 0.04;
    text(ax, 100 .* labelEpsilon, labelRatio, ...
        sprintf("{\\itkh} = %g", khValues(kk)), ...
        "FontSize", cfg.annotationFontSize, "Color", cfg.curveColour, ...
        "HorizontalAlignment", "left", "VerticalAlignment", "bottom", ...
        "Interpreter", "tex");
end
xlabel(ax, "$\epsilon$ (\%)", ...
    "Interpreter", "latex", "FontSize", cfg.labelFontSize);
ylabel(ax, "$h_{\mathrm{est}}/h$", ...
    "Interpreter", "latex", "FontSize", cfg.labelFontSize);
add_panel_label(ax, "(b)", cfg.panelLabelFontSize, [0.97 0.97], "top", "right");
format_axes(ax, cfg);
end

function add_panel_label(ax, labelText, fontSize, position, ...
        verticalAlignment, horizontalAlignment)
if nargin < 6 || strlength(string(horizontalAlignment)) == 0
    horizontalAlignment = "left";
end
text(ax, position(1), position(2), labelText, "Units", "normalized", ...
    "HorizontalAlignment", horizontalAlignment, ...
    "VerticalAlignment", verticalAlignment, "FontSize", fontSize, ...
    "FontWeight", "bold", "Interpreter", "none", "Clipping", "off");
end

function position_curve_legend(ax, legendHandle)
drawnow;
oldAxesUnits = ax.Units;
oldLegendUnits = legendHandle.Units;
ax.Units = "normalized";
legendHandle.Units = "normalized";
axPosition = ax.Position;
legendPosition = legendHandle.Position;
legendPosition(1) = axPosition(1) + axPosition(3) - ...
    legendPosition(3) - 0.04 .* axPosition(3);
legendPosition(2) = axPosition(2) + axPosition(4) - ...
    legendPosition(4) - 0.05 .* axPosition(4);
legendHandle.Position = legendPosition;
ax.Units = oldAxesUnits;
legendHandle.Units = oldLegendUnits;
end

function draw_xy_iqr(ax, x, y, xLow, xHigh, yLow, yHigh, ...
        errorColour, limits)
if isempty(x)
    return
end
capHalfWidth = 0.012 .* diff(xlim(ax));
capHalfHeight = 0.012 .* diff(limits);
for ii = 1:numel(x)
    xEnds = [x(ii)-xLow(ii), x(ii)+xHigh(ii)];
    yEnds = [y(ii)-yLow(ii), y(ii)+yHigh(ii)];
    plot(ax, xEnds, [y(ii) y(ii)], "Color", errorColour, ...
        "LineWidth", 0.7, "HandleVisibility", "off");
    plot(ax, [xEnds(1) xEnds(1)], ...
        [y(ii)-capHalfHeight y(ii)+capHalfHeight], ...
        "Color", errorColour, "LineWidth", 0.7, "HandleVisibility", "off");
    plot(ax, [xEnds(2) xEnds(2)], ...
        [y(ii)-capHalfHeight y(ii)+capHalfHeight], ...
        "Color", errorColour, "LineWidth", 0.7, "HandleVisibility", "off");
    plot(ax, [x(ii) x(ii)], yEnds, "Color", errorColour, ...
        "LineWidth", 0.7, "HandleVisibility", "off");
    xCap = [x(ii)-capHalfWidth x(ii)+capHalfWidth];
    plot(ax, xCap, [yEnds(1) yEnds(1)], "Color", errorColour, ...
        "LineWidth", 0.7, "HandleVisibility", "off");
    plot(ax, xCap, [yEnds(2) yEnds(2)], "Color", errorColour, ...
        "LineWidth", 0.7, "HandleVisibility", "off");
end
end

function limits = padded_limits(lowerValues, upperValues)
lowerValues = lowerValues(isfinite(lowerValues));
upperValues = upperValues(isfinite(upperValues));
if isempty(lowerValues) || isempty(upperValues)
    limits = [0 1];
    return
end
lower = min(lowerValues);
upper = max(upperValues);
if upper <= lower
    upper = lower + 1;
end
padding = 0.07 .* (upper-lower);
limits = [lower-padding upper+padding];
end

function format_axes(ax, cfg)
axis(ax, "square");
set(ax, "FontSize", cfg.axesFontSize, "LineWidth", 0.8, ...
    "TickDir", "out", "TickLabelInterpreter", "tex", "Layer", "top");
end

function map = local_coolwarm(m)
if nargin < 1 || isempty(m)
    m = 256;
end
if m <= 0
    map = zeros(0, 3);
    return
end
anchors = [59 76 192; 84 112 222; 129 164 251; 180 205 251; ...
    221 221 221; 241 184 156; 229 112 88; 203 62 56; 180 4 38] ./ 255;
map = interp1(linspace(0, 1, size(anchors, 1)), anchors, ...
    linspace(0, 1, m), "linear");
map = max(0, min(1, map));
end
