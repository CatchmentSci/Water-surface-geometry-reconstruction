function outputs = generate_figure_08(archiveRoot, outputFolder)
%GENERATE_FIGURE_08 Reproduce the initial-planar-velocity figure.
%   The bundled compact summary uses velocityOutTracked.start.
%   u_streamwise_mps, obtained by projecting tracked paths onto the initial
%   planar water surface. ARCHIVEROOT is retained for API compatibility;
%   the bundled, version-controlled summary always takes precedence.

arguments
    archiveRoot (1, 1) string = "" %#ok<INUSA>
    outputFolder (1, 1) string = fullfile(fileparts(mfilename('fullpath')), "output")
end

outputFolder = string(outputFolder);
summaryFile = fullfile(fileparts(mfilename("fullpath")), "data", ...
    "wse_autocorrelation_velocity_method_sensitivity_summary.mat");

if ~isfile(summaryFile)
    error("Figure08:MissingSummary", ...
        "Required bundled initial-velocity summary was not found:\n%s", summaryFile);
end

source = load(summaryFile, ...
    "autocorrCaseSummaryTable", ...
    "autocorrSensitivitySummaryTable", ...
    "autocorrPanelIncluded");
requiredVariables = ["autocorrCaseSummaryTable", ...
    "autocorrSensitivitySummaryTable", "autocorrPanelIncluded"];
missingVariables = setdiff(requiredVariables, string(fieldnames(source)));
if ~isempty(missingVariables)
    error("Figure08:InvalidSummary", ...
        "Summary file is missing variable(s): %s", ...
        strjoin(missingVariables, ", "));
end

caseSummary = source.autocorrCaseSummaryTable;
sensitivitySummary = source.autocorrSensitivitySummaryTable;
includeMask = logical(source.autocorrPanelIncluded(:));
if ~istable(caseSummary) || ~istable(sensitivitySummary)
    error("Figure08:InvalidSummary", ...
        "The case and sensitivity summaries must be MATLAB tables.");
end
if height(caseSummary) ~= 13 || height(sensitivitySummary) ~= 13 || ...
        numel(includeMask) ~= 13
    error("Figure08:UnexpectedCases", ...
        "Expected the complete 13-case R1-R13 summary.");
end
expectedLabels = "R" + string((1:13).');
if ~isequal(string(caseSummary.caseLabel), expectedLabels) || ...
        ~isequal(string(sensitivitySummary.caseLabel), expectedLabels)
    error("Figure08:UnexpectedCases", ...
        "Summary tables are not the expected ordered R1-R13 set.");
end
if nnz(includeMask) ~= 11 || any(includeMask([12 13]))
    error("Figure08:UnexpectedExclusions", ...
        "Expected 11 plotted cases with R12 and R13 excluded.");
end

requiredCaseFields = ["discharge_Q_m3ps", "deep_nPairs", ...
    "deep_xMedian", "deep_xErrLow", "deep_xErrHigh", ...
    "deep_yMedian", "deep_yErrLow", "deep_yErrHigh"];
requiredSensitivityFields = ["khMedian", "khErrLow", "khErrHigh", ...
    "constant_nPairs", "constant_differenceMedian_pct", ...
    "constant_differenceErrLow_pct", "constant_differenceErrHigh_pct", ...
    "linear_nPairs", "linear_differenceMedian_pct", ...
    "linear_differenceErrLow_pct", "linear_differenceErrHigh_pct", ...
    "power_nPairs", "power_differenceMedian_pct", ...
    "power_differenceErrLow_pct", "power_differenceErrHigh_pct"];
assert_fields(caseSummary, requiredCaseFields, "case summary");
assert_fields(sensitivitySummary, requiredSensitivityFields, ...
    "sensitivity summary");

figCfg = struct;
figCfg.figureWidth_in = 5.5;
figCfg.figureHeight_in = 3.15;
figCfg.figureResolution_dpi = 600;
figCfg.axesFontSize = 9;
figCfg.labelFontSize = 10;
figCfg.panelLabelFontSize = 10;
figCfg.colorbarFontSize = 9;
figCfg.markerArea = 30;
figCfg.errorColour = [0.62 0.62 0.62];

methods = [ ...
    struct("key", "constant", "label", "Constant profile", "marker", "o"), ...
    struct("key", "linear", "label", "Linear profile", "marker", "s"), ...
    struct("key", "power", "label", "Power profile", "marker", "d")];

q = double(caseSummary.discharge_Q_m3ps);
finiteQ = q(isfinite(q) & includeMask);
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
outputBase = fullfile(outputFolder, ...
    "wse_autocorrelation_velocity_validation_and_sensitivity");

fig = figure("Color", "w", "Units", "inches", ...
    "Position", [1 1 figCfg.figureWidth_in figCfg.figureHeight_in], ...
    "PaperPositionMode", "auto");
cleanup = onCleanup(@() close(fig));
layout = tiledlayout(fig, 1, 2, ...
    "TileSpacing", "compact", "Padding", "compact");
colormap(fig, local_coolwarm(256));

ax1 = nexttile(layout, 1);
plot_deep_water_validation(ax1, caseSummary, q, includeMask, ...
    qLimits, figCfg);

ax2 = nexttile(layout, 2);
legendHandle = plot_method_sensitivity(ax2, sensitivitySummary, ...
    methods, q, includeMask, qLimits, figCfg);

cb = colorbar(ax2, "southoutside");
cb.Layout.Tile = "south";
cb.Label.Interpreter = "latex";
cb.Label.String = "$Q$ ($\mathrm{m}^{3}\,\mathrm{s}^{-1}$)";
cb.Label.FontSize = figCfg.labelFontSize;
cb.FontSize = figCfg.colorbarFontSize;
cb.TickLabelInterpreter = "tex";
cb.TickDirection = "out";
position_method_legend(ax2, legendHandle);

pngFile = outputBase + ".png";
pdfFile = outputBase + ".pdf";
exportgraphics(fig, pngFile, "Resolution", figCfg.figureResolution_dpi);
exportgraphics(fig, pdfFile, "ContentType", "vector");

outputs = struct;
outputs.summaryFile = summaryFile;
outputs.pngFile = pngFile;
outputs.pdfFile = pdfFile;
outputs.includedCases = string(caseSummary.caseLabel(includeMask));
outputs.excludedCases = string(caseSummary.caseLabel(~includeMask));
end

function assert_fields(T, requiredFields, description)
missingFields = setdiff(requiredFields, string(T.Properties.VariableNames));
if ~isempty(missingFields)
    error("Figure08:InvalidSummary", "%s is missing field(s): %s", ...
        description, strjoin(missingFields, ", "));
end
end

function plot_deep_water_validation(ax, T, q, includeMask, qLimits, cfg)
x = double(T.deep_xMedian);
y = double(T.deep_yMedian);
xLow = double(T.deep_xErrLow);
xHigh = double(T.deep_xErrHigh);
yLow = double(T.deep_yErrLow);
yHigh = double(T.deep_yErrHigh);
nPairs = double(T.deep_nPairs);
includeMask = logical(includeMask(:));

valid = includeMask & nPairs > 0 & isfinite(x) & isfinite(y) & ...
    isfinite(xLow) & isfinite(xHigh) & isfinite(yLow) & ...
    isfinite(yHigh) & isfinite(q);
x = x(valid);
y = y(valid);
xLow = xLow(valid);
xHigh = xHigh(valid);
yLow = yLow(valid);
yHigh = yHigh(valid);
q = q(valid);

limits = padded_limits([x-xLow; y-yLow], [x+xHigh; y+yHigh], true);
hold(ax, "on");
grid(ax, "off");
box(ax, "on");
xlim(ax, limits);
ylim(ax, limits);
plot(ax, limits, limits, "--", "Color", [0.20 0.20 0.20], ...
    "LineWidth", 0.8, "HandleVisibility", "off");
draw_xy_iqr(ax, x, y, xLow, xHigh, yLow, yHigh, ...
    cfg.errorColour, limits);
scatter(ax, x, y, cfg.markerArea, q, "filled", ...
    "MarkerEdgeColor", [0.15 0.15 0.15], "LineWidth", 0.5);
clim(ax, qLimits);
xlabel(ax, "$U_s$ ($\mathrm{m}\,\mathrm{s}^{-1}$)", ...
    "Interpreter", "latex", "FontSize", cfg.labelFontSize);
ylabel(ax, "$U_{\mathrm{deep}}$ ($\mathrm{m}\,\mathrm{s}^{-1}$)", ...
    "Interpreter", "latex", "FontSize", cfg.labelFontSize);
add_panel_label(ax, "(a)", cfg.panelLabelFontSize, [0.03 0.97], "top");
format_axes(ax, cfg);
end

function legendHandle = plot_method_sensitivity( ...
        ax, T, methods, q, includeMask, qLimits, cfg)
x = double(T.khMedian);
xLow = double(T.khErrLow);
xHigh = double(T.khErrHigh);
includeMask = logical(includeMask(:));

allLowerY = 0;
allUpperY = 0;
for mm = 1:numel(methods)
    key = methods(mm).key;
    y = double(T.(key + "_differenceMedian_pct"));
    yLow = double(T.(key + "_differenceErrLow_pct"));
    yHigh = double(T.(key + "_differenceErrHigh_pct"));
    validY = includeMask & isfinite(y) & isfinite(yLow) & isfinite(yHigh);
    allLowerY = [allLowerY; y(validY)-yLow(validY)]; %#ok<AGROW>
    allUpperY = [allUpperY; y(validY)+yHigh(validY)]; %#ok<AGROW>
end

validX = includeMask & isfinite(x) & x > 0 & isfinite(xLow) & ...
    isfinite(xHigh) & (x-xLow) > 0 & isfinite(q);
xLimits = padded_limits(x(validX)-xLow(validX), ...
    x(validX)+xHigh(validX), false);
yLimits = padded_limits(allLowerY, allUpperY, false);

hold(ax, "on");
grid(ax, "off");
box(ax, "on");
xlim(ax, xLimits);
ylim(ax, yLimits);
legendHandles = gobjects(1, numel(methods));
for mm = 1:numel(methods)
    key = methods(mm).key;
    y = double(T.(key + "_differenceMedian_pct"));
    yLow = double(T.(key + "_differenceErrLow_pct"));
    yHigh = double(T.(key + "_differenceErrHigh_pct"));
    nPairs = double(T.(key + "_nPairs"));
    valid = validX & nPairs > 0 & isfinite(y) & ...
        isfinite(yLow) & isfinite(yHigh);
    draw_xy_iqr(ax, x(valid), y(valid), xLow(valid), xHigh(valid), ...
        yLow(valid), yHigh(valid), cfg.errorColour, yLimits);
    scatter(ax, x(valid), y(valid), cfg.markerArea, q(valid), ...
        methods(mm).marker, "filled", ...
        "MarkerEdgeColor", [0.15 0.15 0.15], "LineWidth", 0.5);
    legendHandles(mm) = scatter(ax, nan, nan, cfg.markerArea, ...
        [0.72 0.72 0.72], methods(mm).marker, "filled", ...
        "MarkerEdgeColor", [0.15 0.15 0.15], "LineWidth", 0.5);
end
clim(ax, qLimits);
legendHandle = legend(ax, legendHandles, {methods.label}, ...
    "Location", "none", "FontSize", 7.5, "Box", "off", ...
    "Interpreter", "none");
xlabel(ax, "$kh=2\pi h/\lambda$", ...
    "Interpreter", "latex", "FontSize", cfg.labelFontSize);
ylabel(ax, ...
    "$100\,(U_{\mathrm{method}}-U_{\mathrm{deep}})/U_{\mathrm{deep}}$ (\%)", ...
    "Interpreter", "latex", "FontSize", cfg.labelFontSize);
add_panel_label(ax, "(b)", cfg.panelLabelFontSize, ...
    [0.97 0.97], "top", "right");
format_axes(ax, cfg);
end

function position_method_legend(ax, legendHandle)
drawnow;
oldAxesUnits = ax.Units;
oldLegendUnits = legendHandle.Units;
ax.Units = "normalized";
legendHandle.Units = "normalized";
axPosition = ax.Position;
legendPosition = legendHandle.Position;
legendPosition(1) = min(axPosition(1) + 0.50 .* axPosition(3), ...
    axPosition(1) + axPosition(3) - legendPosition(3) - 0.03 .* axPosition(3));
legendPosition(2) = axPosition(2) + axPosition(4) - ...
    legendPosition(4) - 0.12 .* axPosition(4);
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

function limits = padded_limits(lowerValues, upperValues, equalAxes)
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
if equalAxes && limits(1) < 0 && lower >= 0
    limits(1) = 0;
end
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
