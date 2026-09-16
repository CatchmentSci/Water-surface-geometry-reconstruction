function outputs = generate_figure_09(archiveRoot, outputFolder)
%GENERATE_FIGURE_09 Reproduce the initial-planar-velocity depth figure.
%   The bundled depth summary was re-inverted using
%   velocityOutTracked.start.u_streamwise_mps. ARCHIVEROOT is retained for
%   API compatibility; the bundled, version-controlled summary always
%   takes precedence.

arguments
    archiveRoot (1, 1) string = "" %#ok<INUSA>
    outputFolder (1, 1) string = fullfile(fileparts(mfilename('fullpath')), "output")
end

outputFolder = string(outputFolder);
summaryFile = fullfile(fileparts(mfilename("fullpath")), "data", ...
    "wse_autocorrelation_uniform_linear_power_depth_summary.csv");
if ~isfile(summaryFile)
    error("Figure09:MissingSummary", ...
        "Required bundled initial-velocity summary was not found:\n%s", summaryFile);
end

T = readtable(summaryFile, "VariableNamingRule", "preserve");
if height(T) ~= 13
    error("Figure09:UnexpectedCases", ...
        "Expected the complete 13-case R1-R13 summary.");
end
expectedLabels = "R" + string((1:13).');
if ~isequal(string(T.caseLabel), expectedLabels)
    error("Figure09:UnexpectedCases", ...
        "Summary table is not the expected ordered R1-R13 set.");
end

methods = [ ...
    struct("key", "uniform", "panelLabel", "(a)"), ...
    struct("key", "linear", "panelLabel", "(b)"), ...
    struct("key", "power", "panelLabel", "(c)")];
branches = [struct("key", "shallow"), struct("key", "deep")];
requiredFields = ["caseLabel", "discharge_Q_m3ps", "includedInFigure"];
for mm = 1:numel(methods)
    for bb = 1:numel(branches)
        prefix = methods(mm).key + "_" + branches(bb).key;
        requiredFields = [requiredFields, prefix + [ ...
            "_nPairs", "_xMedian", "_xErrLow", "_xErrHigh", ...
            "_yMedian", "_yErrLow", "_yErrHigh"]]; %#ok<AGROW>
    end
end
missingFields = setdiff(requiredFields, string(T.Properties.VariableNames));
if ~isempty(missingFields)
    error("Figure09:InvalidSummary", ...
        "Summary table is missing field(s): %s", ...
        strjoin(missingFields, ", "));
end

included = logical(T.includedInFigure(:));
if nnz(included) ~= 11 || ~all(included(1:11)) || any(included(12:13))
    error("Figure09:UnexpectedExclusions", ...
        "Expected 11 plotted cases with R12 and R13 excluded.");
end

cfg = struct;
cfg.figureWidth_in = 6.6;
cfg.figureHeight_in = 3.00;
cfg.figureResolution_dpi = 600;
cfg.axesFontSize = 8.5;
cfg.labelFontSize = 9.5;
cfg.panelLabelFontSize = 9;
cfg.colorbarFontSize = 8.5;
cfg.markerArea = 28;
cfg.errorColour = [0.62 0.62 0.62];

q = double(T.discharge_Q_m3ps);
finiteQ = q(isfinite(q) & included);
if isempty(finiteQ)
    qLimits = [0 1];
else
    qLimits = [min(finiteQ), max(finiteQ)];
    if qLimits(1) == qLimits(2)
        qLimits = qLimits + [-0.5 0.5];
    end
end

% The paper uses identical limits for shallow- and deep-solution figures.
commonAxisLimits = calculate_common_axis_limits( ...
    T, methods, branches, included);

if ~isfolder(outputFolder)
    mkdir(outputFolder);
end
outputBase = fullfile(outputFolder, ...
    "Figure9");

fig = make_depth_figure(T, methods, included, q, qLimits, ...
    commonAxisLimits, cfg);
cleanup = onCleanup(@() close(fig));
pngFile = outputBase + ".png";
pdfFile = outputBase + ".pdf";
exportgraphics(fig, pngFile, "Resolution", cfg.figureResolution_dpi);
exportgraphics(fig, pdfFile, "ContentType", "vector");

outputs = struct;
outputs.summaryFile = summaryFile;
outputs.pngFile = pngFile;
outputs.pdfFile = pdfFile;
outputs.includedCases = string(T.caseLabel(included));
outputs.excludedCases = string(T.caseLabel(~included));
outputs.axisLimits_m = commonAxisLimits;
end

function limits = calculate_common_axis_limits(T, methods, branches, included)
lowerValues = [];
upperValues = [];
for mm = 1:numel(methods)
    for bb = 1:numel(branches)
        prefix = methods(mm).key + "_" + branches(bb).key;
        valid = included & T.(prefix + "_nPairs") > 0;
        x = double(T.(prefix + "_xMedian"));
        y = double(T.(prefix + "_yMedian"));
        xLow = double(T.(prefix + "_xErrLow"));
        xHigh = double(T.(prefix + "_xErrHigh"));
        yLow = double(T.(prefix + "_yErrLow"));
        yHigh = double(T.(prefix + "_yErrHigh"));
        lowerValues = [lowerValues; x(valid)-xLow(valid); ...
            y(valid)-yLow(valid)]; %#ok<AGROW>
        upperValues = [upperValues; x(valid)+xHigh(valid); ...
            y(valid)+yHigh(valid)]; %#ok<AGROW>
    end
end
lowerValues = lowerValues(isfinite(lowerValues));
upperValues = upperValues(isfinite(upperValues));
if isempty(lowerValues) || isempty(upperValues)
    limits = [0 1];
    return
end
lower = min(lowerValues);
upper = max(upperValues);
span = upper - lower;
if span <= 0
    span = 1;
end
limits = [max(0, lower-0.07*span), upper+0.07*span];
end

function fig = make_depth_figure(T, methods, included, q, qLimits, limits, cfg)
fig = figure("Color", "w", "Name", "Deep-depth solutions", ...
    "NumberTitle", "off", "Units", "inches", ...
    "Position", [1 1 cfg.figureWidth_in cfg.figureHeight_in], ...
    "PaperPositionMode", "auto");
layout = tiledlayout(fig, 1, 3, ...
    "TileSpacing", "compact", "Padding", "compact");
colormap(fig, local_coolwarm(256));
ax = gobjects(1, numel(methods));
for mm = 1:numel(methods)
    ax(mm) = nexttile(layout, mm);
    plot_depth_panel(ax(mm), T, methods(mm), included, ...
        q, qLimits, limits, cfg);
    if mm > 1
        ax(mm).YTickLabel = [];
    end
end
xlabel(ax(2), "$h$ ($\mathrm{m}$)", ...
    "Interpreter", "latex", "FontSize", cfg.labelFontSize);
ylabel(ax(1), "$h_{\mathrm{est}}$ ($\mathrm{m}$)", ...
    "Interpreter", "latex", "FontSize", cfg.labelFontSize);
cb = colorbar(ax(end), "southoutside");
cb.Layout.Tile = "south";
cb.Label.Interpreter = "latex";
cb.Label.String = "$Q$ ($\mathrm{m}^{3}\,\mathrm{s}^{-1}$)";
cb.Label.FontSize = cfg.labelFontSize;
cb.FontSize = cfg.colorbarFontSize;
cb.TickLabelInterpreter = "tex";
cb.TickDirection = "out";
end

function plot_depth_panel(ax, T, method, included, q, qLimits, limits, cfg)
prefix = method.key + "_deep";
nPairs = double(T.(prefix + "_nPairs"));
x = double(T.(prefix + "_xMedian"));
y = double(T.(prefix + "_yMedian"));
xLow = double(T.(prefix + "_xErrLow"));
xHigh = double(T.(prefix + "_xErrHigh"));
yLow = double(T.(prefix + "_yErrLow"));
yHigh = double(T.(prefix + "_yErrHigh"));
valid = included & nPairs > 0 & isfinite(x) & isfinite(y) & ...
    isfinite(xLow) & isfinite(xHigh) & isfinite(yLow) & ...
    isfinite(yHigh) & isfinite(q);
x = x(valid);
y = y(valid);
xLow = xLow(valid);
xHigh = xHigh(valid);
yLow = yLow(valid);
yHigh = yHigh(valid);
q = q(valid);

hold(ax, "on");
grid(ax, "off");
box(ax, "on");
xlim(ax, limits);
ylim(ax, limits);
plot(ax, limits, limits, "--", "Color", [0.20 0.20 0.20], ...
    "LineWidth", 0.8, "HandleVisibility", "off");
draw_xy_iqr(ax, x, y, xLow, xHigh, yLow, yHigh, cfg.errorColour);
scatter(ax, x, y, cfg.markerArea, q, "filled", ...
    "MarkerEdgeColor", [0.15 0.15 0.15], "LineWidth", 0.5);
clim(ax, qLimits);
text(ax, 0.03, 0.97, method.panelLabel, "Units", "normalized", ...
    "HorizontalAlignment", "left", "VerticalAlignment", "top", ...
    "FontSize", cfg.panelLabelFontSize, "FontWeight", "bold", ...
    "Interpreter", "none", "Clipping", "off");
axis(ax, "square");
set(ax, "FontSize", cfg.axesFontSize, "LineWidth", 0.8, ...
    "TickDir", "out", "TickLabelInterpreter", "tex", "Layer", "top");
end

function draw_xy_iqr(ax, x, y, xLow, xHigh, yLow, yHigh, colour)
if isempty(x)
    return
end
capHalfWidth = 0.012 .* diff(xlim(ax));
capHalfHeight = 0.012 .* diff(ylim(ax));
for ii = 1:numel(x)
    xEnds = [x(ii)-xLow(ii), x(ii)+xHigh(ii)];
    yEnds = [y(ii)-yLow(ii), y(ii)+yHigh(ii)];
    plot(ax, xEnds, [y(ii) y(ii)], "Color", colour, ...
        "LineWidth", 0.7, "HandleVisibility", "off");
    plot(ax, [xEnds(1) xEnds(1)], ...
        [y(ii)-capHalfHeight y(ii)+capHalfHeight], ...
        "Color", colour, "LineWidth", 0.7, "HandleVisibility", "off");
    plot(ax, [xEnds(2) xEnds(2)], ...
        [y(ii)-capHalfHeight y(ii)+capHalfHeight], ...
        "Color", colour, "LineWidth", 0.7, "HandleVisibility", "off");
    plot(ax, [x(ii) x(ii)], yEnds, "Color", colour, ...
        "LineWidth", 0.7, "HandleVisibility", "off");
    plot(ax, [x(ii)-capHalfWidth x(ii)+capHalfWidth], ...
        [yEnds(1) yEnds(1)], "Color", colour, ...
        "LineWidth", 0.7, "HandleVisibility", "off");
    plot(ax, [x(ii)-capHalfWidth x(ii)+capHalfWidth], ...
        [yEnds(2) yEnds(2)], "Color", colour, ...
        "LineWidth", 0.7, "HandleVisibility", "off");
end
end

function map = local_coolwarm(m)
if nargin < 1 || isempty(m)
    m = 256;
end
anchors = [59 76 192; 84 112 222; 129 164 251; 180 205 251; ...
    221 221 221; 241 184 156; 229 112 88; 203 62 56; 180 4 38] ./ 255;
map = interp1(linspace(0, 1, size(anchors, 1)), anchors, ...
    linspace(0, 1, m), "linear");
map = max(0, min(1, map));
end
