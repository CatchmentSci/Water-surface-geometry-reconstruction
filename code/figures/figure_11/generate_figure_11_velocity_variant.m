function outputs = generate_figure_11_velocity_variant(archiveRoot, batchMatFile, outputFolder, options)
%GENERATE_FIGURE_11_VELOCITY_VARIANT Figure 11 with a chosen tracked-velocity source.
%   OUTPUTS = GENERATE_FIGURE_11_VELOCITY_VARIANT(ARCHIVEROOT, BATCHMATFILE,
%   OUTPUTFOLDER) rebuilds the observed points of Figure 11(a) from a
%   transect-level depth inversion driven by the tracked surface velocity
%   projected onto the ACCEPTED water-surface map, rather than the velocity
%   deposited with the archived depth summary, and writes the two-panel
%   figure, a per-case CSV comparing the velocity sources, and a
%   per-transect CSV.
%
%   Inputs
%     ARCHIVEROOT   root of the extracted Zenodo archive
%     BATCHMATFILE  optional path to
%                   checkpoint_batch_initial_vs_accepted_velocity_summary.mat
%                   produced by batch_first10m_initial_vs_accepted_velocity_scatterplots.m;
%                   defaults to the originating workstation path shown in
%                   the arguments block below
%     OUTPUTFOLDER  destination (default: ./output beside this file)
%
%   Per transect, with k = 2*pi/lambda from the accepted-map autocorrelation
%   wavelength of the Figure 8 summary and U the chosen tracked velocity,
%   the constant-profile inversion (Eq. 10a) has the single admissible root
%
%       h_est = atanh(U^2 k / g) / k,      valid only where U^2 k / g < 1 .
%
%   Case values are paired medians of surveyed depth h and h_est over the
%   transects that are commonValid in the Figure 8 summary and admissible
%   under the chosen velocity, exactly as for Figure 9(a). Panel (a) plots
%   h_est/h against kh with the predicted amplification curves; panel (b)
%   is unchanged from generate_figure_11.
%
%   Name-value options
%     "VelocitySource"   "accepted" (default) | "initial" | "archived"
%                        "archived" reproduces the deposited Figure 9/11
%                        values from the archive CSVs without re-inversion
%                        and is the reference the other two are compared to.
%     "EpsilonValues"    panel (a) curves (default [0.08 0.12 0.16])
%     "KhValues"         panel (b) curves (default [1 3 8])
%     "ReferenceEpsilon" panel (b) marker (default 0.08)
%     "MaxDepth_m"       discard transect inversions deeper than this
%                        (default Inf; set e.g. 10 to mimic a "reasonable
%                        range" filter if the archived analysis used one)
%
%   The function also reports which of the two batch velocities the
%   archived per-transect velocity matches, so the provenance of the
%   deposited figure is checked rather than assumed.

arguments
    archiveRoot (1, 1) string
    batchMatFile (1, 1) string = ...
        "D:\OneDrive - Newcastle University\Documents - WSE Project\General\Dart\Videos\Inputs\batch_first10m_initial_vs_accepted_velocity_outputs\checkpoint_batch_initial_vs_accepted_velocity_summary.mat"
    outputFolder (1, 1) string = fullfile(fileparts(mfilename('fullpath')), "output")
    options.VelocitySource (1, 1) string {mustBeMember(options.VelocitySource, ["accepted", "initial", "archived"])} = "accepted"
    options.EpsilonValues (1, :) double {mustBePositive, mustBeLessThan(options.EpsilonValues, 1)} = [0.08 0.12 0.16]
    options.KhValues (1, :) double {mustBePositive} = [1 3 8]
    options.ReferenceEpsilon (1, 1) double {mustBeNonnegative, mustBeLessThan(options.ReferenceEpsilon, 1)} = 0.08
    options.MaxDepth_m (1, 1) double {mustBePositive} = Inf
end

g = 9.81;
archiveRoot = string(archiveRoot);
outputRoot = fullfile(archiveRoot, "videos", "outputs");
depthFile = fullfile(outputRoot, "wse_autocorrelation_uniform_linear_power_depth_summary.csv");
velocityFile = fullfile(outputRoot, "wse_autocorrelation_velocity_method_sensitivity_summary.mat");
requiredFiles = [depthFile, velocityFile, string(batchMatFile)];
for ii = 1:numel(requiredFiles)
    if ~isfile(requiredFiles(ii))
        error("Figure11Variant:MissingInput", ...
            "Required file was not found:\n%s", requiredFiles(ii));
    end
end

%% ---- archived case-level summary (reference)
depth = readtable(depthFile, "VariableNamingRule", "preserve");
expectedLabels = "R" + string((1:13).');
if height(depth) ~= 13 || ~isequal(string(depth.caseLabel), expectedLabels)
    error("Figure11Variant:UnexpectedCases", "Depth summary is not the ordered R1-R13 set.");
end
included = logical(depth.includedInFigure(:));
q = double(depth.discharge_Q_m3ps);

%% ---- Figure 8 per-transect tables (accepted-map autocorrelation wavelength)
fig8 = load(velocityFile, "autocorrCaseTransectTables", "autocorrPanelIncluded");
transects = fig8.autocorrCaseTransectTables;
if numel(transects) ~= 13 || ~isequal(logical(fig8.autocorrPanelIncluded(:)), included)
    error("Figure11Variant:InclusionMismatch", ...
        "Figure 8 transect tables do not carry the expected 13 cases and R1-R11 mask.");
end

%% ---- batch initial/accepted velocities (per transect)
batch = load(batchMatFile, "caseOutputs", "caseSummaryTable");
caseOutputs = batch.caseOutputs;
batchSummary = batch.caseSummaryTable;
if numel(caseOutputs) ~= 13 || ~isequal(string(batchSummary.caseLabel), expectedLabels)
    error("Figure11Variant:UnexpectedCases", "Batch MAT is not the ordered R1-R13 set.");
end

%% ---- transect-level assembly
perTransect = table();
caseRows = repmat(empty_case_row(), 13, 1);
for ii = 1:13
    T8 = transects{ii};
    assert_fields(T8, ["caseLabel", "acceptedAutocorrWavelength_m", "commonValid"], ...
        sprintf("Figure 8 transect table %d", ii));
    if any(string(T8.caseLabel) ~= expectedLabels(ii))
        error("Figure11Variant:CaseAlignment", "Transect table %d is not %s.", ii, expectedLabels(ii));
    end

    A = caseOutputs{ii}.transectAnalysis;
    if ~isstruct(A) || isempty(A)
        error("Figure11Variant:MissingTransects", "Batch transectAnalysis is empty for %s.", expectedLabels(ii));
    end
    profileIndex = [A.profileIndex].';
    rowCoord = [A.rowCoord].';
    uInitial = [A.initialVelocityTrackedMedian_mps].';
    uAccepted = [A.acceptedVelocityTrackedMedian_mps].';

    % Align the Figure 8 transect rows to the batch transects.
    keyName = first_field(T8, ["profileIndex", "transectIndex", "rowCoord"]);
    if strlength(keyName) == 0
        if height(T8) ~= numel(profileIndex)
            error("Figure11Variant:TransectAlignment", ...
                "Cannot align Figure 8 transects for %s (no key column, %d vs %d rows).", ...
                expectedLabels(ii), height(T8), numel(profileIndex));
        end
        idx8 = (1:height(T8)).';
    elseif keyName == "rowCoord"
        [tf, idx8] = ismember(round(rowCoord), round(double(T8.rowCoord)));
        if ~all(tf), error("Figure11Variant:TransectAlignment", "rowCoord mismatch for %s.", expectedLabels(ii)); end
    else
        [tf, idx8] = ismember(profileIndex, double(T8.(keyName)));
        if ~all(tf), error("Figure11Variant:TransectAlignment", "%s mismatch for %s.", keyName, expectedLabels(ii)); end
    end
    lambda = double(T8.acceptedAutocorrWavelength_m(idx8));
    commonValid = logical(T8.commonValid(idx8));

    % Surveyed depth per transect from the archived per-transect product.
    baseName = string(caseOutputs{ii}.baseName);
    perTransectFile = fullfile(outputRoot, baseName + ...
        "_first10m_transect_depth_from_observed_velocity_with_deep_water_velocity_pmusic_depths.csv");
    if isfile(perTransectFile)
        P = readtable(perTransectFile, "VariableNamingRule", "preserve");
        assert_fields(P, ["profileIndex", "rowCoord", "crossSectionDepth_m", ...
            "observedVelocityTrackedMedian_mps"], "archived per-transect table");
        [tf, idxP] = ismember(profileIndex, double(P.profileIndex));
        if ~all(tf) || any(abs(double(P.rowCoord(idxP)) - rowCoord) > 0.5)
            error("Figure11Variant:TransectAlignment", ...
                "Archived per-transect rows do not match the batch transects for %s.", expectedLabels(ii));
        end
        hSurvey = double(P.crossSectionDepth_m(idxP));
        uArchivedTransect = double(P.observedVelocityTrackedMedian_mps(idxP));
    elseif included(ii)
        error("Figure11Variant:MissingInput", "Archived per-transect file not found:\n%s", perTransectFile);
    else
        % R12 and R13 are excluded from the figure and have no deposited
        % per-transect product; carry them through as NaN.
        warning("Figure11Variant:MissingInput", ...
            "No archived per-transect file for %s (excluded case); values set to NaN.", expectedLabels(ii));
        hSurvey = nan(size(profileIndex));
        uArchivedTransect = nan(size(profileIndex));
    end

    % Which batch velocity does the archived per-transect velocity match?
    ok = isfinite(uArchivedTransect) & isfinite(uInitial) & isfinite(uAccepted);
    dAcc = max(abs(uArchivedTransect(ok) - uAccepted(ok)));
    dIni = max(abs(uArchivedTransect(ok) - uInitial(ok)));
    if isempty(dAcc), dAcc = NaN; dIni = NaN; end

    % Velocity used by the Figure 8 summary, if it carries one.
    velName = first_field(T8, ["trackedVelocity_mps", "observedVelocity_mps", ...
        "observedVelocityTrackedMedian_mps", "Us_mps", "U_s_mps"]);
    if strlength(velName) > 0
        uFig8 = double(T8.(velName)(idx8));
    else
        uFig8 = nan(size(lambda));
    end

    k = 2 .* pi ./ lambda;
    hEstAccepted = invert_constant_profile(uAccepted, k, g, options.MaxDepth_m);
    hEstInitial = invert_constant_profile(uInitial, k, g, options.MaxDepth_m);
    hEstArchivedTransect = invert_constant_profile(uArchivedTransect, k, g, options.MaxDepth_m);

    n = numel(profileIndex);
    perTransect = [perTransect; table( ...
        repmat(ii, n, 1), repmat(expectedLabels(ii), n, 1), profileIndex, rowCoord, ...
        commonValid, lambda, hSurvey, uInitial, uAccepted, uArchivedTransect, uFig8, ...
        hEstInitial, hEstAccepted, hEstArchivedTransect, ...
        VariableNames=["caseNumber", "caseLabel", "profileIndex", "rowCoord", ...
        "commonValid", "lambda_m", "hSurveyed_m", "U_initial_mps", "U_accepted_mps", ...
        "U_archivedTransect_mps", "U_figure8_mps", ...
        "hEst_initial_m", "hEst_accepted_m", "hEst_archivedTransect_m"])]; %#ok<AGROW>

    r = empty_case_row();
    r.caseNumber = ii;
    r.caseLabel = expectedLabels(ii);
    r.discharge_Q_m3ps = q(ii);
    r.included = included(ii);
    r.archivedMatchesAccepted_maxAbsDiff_mps = dAcc;
    r.archivedMatchesInitial_maxAbsDiff_mps = dIni;
    r.U_initial_median_mps = median(uInitial, "omitnan");
    r.U_accepted_median_mps = median(uAccepted, "omitnan");
    r.U_accepted_over_initial_minus1_pct = 100 .* (r.U_accepted_median_mps ./ r.U_initial_median_mps - 1);
    r.lambda_median_m = median(lambda(commonValid), "omitnan");
    r = add_paired_medians(r, "accepted", hSurvey, hEstAccepted, commonValid, k);
    r = add_paired_medians(r, "initial", hSurvey, hEstInitial, commonValid, k);
    % Reference: deposited Figure 9 case medians (no re-inversion).
    r.archived_nPairs = double(depth.uniform_deep_nPairs(ii));
    r.archived_h_m = double(depth.uniform_deep_xMedian(ii));
    r.archived_hEst_m = double(depth.uniform_deep_yMedian(ii));
    r.archived_kh = 2 .* pi .* r.archived_h_m ./ r.lambda_median_m;
    r.archived_ratio = r.archived_hEst_m ./ r.archived_h_m;
    r.archived_ratioErrLow = r.archived_ratio - double(depth.uniform_deep_yQ25(ii)) ./ r.archived_h_m;
    r.archived_ratioErrHigh = double(depth.uniform_deep_yQ75(ii)) ./ r.archived_h_m - r.archived_ratio;
    r.archived_khErrLow = r.archived_kh - 2 .* pi .* double(depth.uniform_deep_xQ25(ii)) ./ r.lambda_median_m;
    r.archived_khErrHigh = 2 .* pi .* double(depth.uniform_deep_xQ75(ii)) ./ r.lambda_median_m - r.archived_kh;
    caseRows(ii) = r;
end
caseTable = struct2table(caseRows);

% Implied velocity bias for each source.
for src = ["accepted", "initial", "archived"]
    kh = caseTable.(src + "_kh");
    ratio = caseTable.(src + "_ratio");
    caseTable.(src + "_impliedEpsilon") = 1 - sqrt(tanh(ratio .* kh) ./ tanh(kh));
end
caseTable.deltaRatio_acceptedMinusArchived = caseTable.accepted_ratio - caseTable.archived_ratio;
caseTable.deltaRatio_initialMinusArchived = caseTable.initial_ratio - caseTable.archived_ratio;

fprintf("\nArchived per-transect velocity vs batch velocities (max |diff|, m/s):\n");
for ii = 1:13
    fprintf("  %-4s accepted %.4f   initial %.4f\n", expectedLabels(ii), ...
        caseTable.archivedMatchesAccepted_maxAbsDiff_mps(ii), ...
        caseTable.archivedMatchesInitial_maxAbsDiff_mps(ii));
end

%% ---- select the source to plot
src = options.VelocitySource;
plotKh = caseTable.(src + "_kh");
plotRatio = caseTable.(src + "_ratio");
plotKhLow = caseTable.(src + "_khErrLow");
plotKhHigh = caseTable.(src + "_khErrHigh");
plotRatioLow = caseTable.(src + "_ratioErrLow");
plotRatioHigh = caseTable.(src + "_ratioErrHigh");
plotN = caseTable.(src + "_nPairs");
valid = included & plotN > 0 & isfinite(plotKh) & isfinite(plotRatio) & ...
    isfinite(plotKhLow) & isfinite(plotKhHigh) & isfinite(plotRatioLow) & ...
    isfinite(plotRatioHigh) & isfinite(q);

%% ---- figure (Figure 8 house style)
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
    if qLimits(1) == qLimits(2), qLimits = qLimits + [-0.5 0.5]; end
end

if ~isfolder(outputFolder), mkdir(outputFolder); end
outputBase = fullfile(outputFolder, ...
    "wse_autocorrelation_depth_bias_amplification_" + src + "_velocity");

fig = figure("Color", "w", "Units", "inches", ...
    "Position", [1 1 figCfg.figureWidth_in figCfg.figureHeight_in], ...
    "PaperPositionMode", "auto");
cleanup = onCleanup(@() close(fig));
layout = tiledlayout(fig, 1, 2, "TileSpacing", "compact", "Padding", "compact");
colormap(fig, local_coolwarm(256));

ax1 = nexttile(layout, 1);
legendHandle = plot_amplification_panel(ax1, plotKh(valid), plotRatio(valid), ...
    plotKhLow(valid), plotKhHigh(valid), plotRatioLow(valid), plotRatioHigh(valid), ...
    q(valid), qLimits, options.EpsilonValues, figCfg);
ax2 = nexttile(layout, 2);
plot_sensitivity_panel(ax2, options.KhValues, options.ReferenceEpsilon, figCfg);

cb = colorbar(ax1, "southoutside");
cb.Layout.Tile = "south";
cb.Label.Interpreter = "latex";
cb.Label.String = "$Q$ ($\mathrm{m}^{3}\,\mathrm{s}^{-1}$)";
cb.Label.FontSize = figCfg.labelFontSize;
cb.FontSize = figCfg.colorbarFontSize;
cb.TickLabelInterpreter = "tex";
cb.TickDirection = "out";
position_curve_legend(ax1, legendHandle);

exportgraphics(fig, outputBase + ".png", "Resolution", figCfg.figureResolution_dpi);
exportgraphics(fig, outputBase + ".pdf", "ContentType", "vector");
caseCsv = fullfile(outputFolder, "figure_11_velocity_source_comparison_by_case.csv");
transectCsv = fullfile(outputFolder, "figure_11_velocity_source_comparison_by_transect.csv");
writetable(caseTable, caseCsv);
writetable(perTransect, transectCsv);

outputs = struct;
outputs.velocitySource = src;
outputs.pngFile = outputBase + ".png";
outputs.pdfFile = outputBase + ".pdf";
outputs.caseCsv = caseCsv;
outputs.transectCsv = transectCsv;
outputs.caseTable = caseTable;
outputs.perTransect = perTransect;
outputs.plottedCases = expectedLabels(valid);
end

%% ======================================================================
function hEst = invert_constant_profile(U, k, g, maxDepth)
% Single admissible root of tanh(k h) = U^2 k / g (constant profile, Eq. 10a).
x = U.^2 .* k ./ g;
hEst = nan(size(U));
ok = isfinite(x) & isfinite(k) & k > 0 & x > 0 & x < 1;
hEst(ok) = atanh(x(ok)) ./ k(ok);
hEst(hEst > maxDepth) = NaN;
end

function r = add_paired_medians(r, prefix, h, hEst, commonValid, k)
valid = commonValid & isfinite(h) & isfinite(hEst) & isfinite(k);
r.(prefix + "_nPairs") = nnz(valid);
fields = ["_h_m", "_hEst_m", "_kh", "_ratio", "_khErrLow", "_khErrHigh", "_ratioErrLow", "_ratioErrHigh"];
for f = fields, r.(prefix + f) = NaN; end
if ~any(valid), return; end
hv = h(valid); ev = hEst(valid); kv = k(valid);
hMed = median(hv); eMed = median(ev); kMed = median(kv);
r.(prefix + "_h_m") = hMed;
r.(prefix + "_hEst_m") = eMed;
r.(prefix + "_kh") = kMed .* hMed;
r.(prefix + "_ratio") = eMed ./ hMed;
r.(prefix + "_khErrLow") = kMed .* (hMed - prctile(hv, 25));
r.(prefix + "_khErrHigh") = kMed .* (prctile(hv, 75) - hMed);
r.(prefix + "_ratioErrLow") = (eMed - prctile(ev, 25)) ./ hMed;
r.(prefix + "_ratioErrHigh") = (prctile(ev, 75) - eMed) ./ hMed;
end

function r = empty_case_row()
r = struct("caseNumber", NaN, "caseLabel", "", "discharge_Q_m3ps", NaN, "included", false, ...
    "archivedMatchesAccepted_maxAbsDiff_mps", NaN, "archivedMatchesInitial_maxAbsDiff_mps", NaN, ...
    "U_initial_median_mps", NaN, "U_accepted_median_mps", NaN, ...
    "U_accepted_over_initial_minus1_pct", NaN, "lambda_median_m", NaN);
for prefix = ["accepted", "initial"]
    r.(prefix + "_nPairs") = 0;
    for f = ["_h_m", "_hEst_m", "_kh", "_ratio", "_khErrLow", "_khErrHigh", "_ratioErrLow", "_ratioErrHigh"]
        r.(prefix + f) = NaN;
    end
end
for f = ["archived_nPairs", "archived_h_m", "archived_hEst_m", "archived_kh", "archived_ratio", ...
        "archived_ratioErrLow", "archived_ratioErrHigh", "archived_khErrLow", "archived_khErrHigh"]
    r.(f) = NaN;
end
end

function name = first_field(T, candidates)
name = "";
vars = string(T.Properties.VariableNames);
for c = candidates
    if any(vars == c), name = c; return; end
end
end

function assert_fields(T, requiredFields, description)
missingFields = setdiff(requiredFields, string(T.Properties.VariableNames));
if ~isempty(missingFields)
    error("Figure11Variant:InvalidInput", "%s is missing field(s): %s", ...
        description, strjoin(missingFields, ", "));
end
end

function r = predicted_ratio(kh, epsilon)
r = atanh((1 - epsilon).^2 .* tanh(kh)) ./ kh;
end

function legendHandle = plot_amplification_panel(ax, kh, ratio, ...
        khErrLow, khErrHigh, ratioErrLow, ratioErrHigh, q, qLimits, epsilonValues, cfg)
xLimits = padded_limits(kh - khErrLow, kh + khErrHigh); xLimits(1) = 0;
yLimits = padded_limits(ratio - ratioErrLow, ratio + ratioErrHigh); yLimits(1) = 0;
yLimits(2) = max(yLimits(2), 1.0);
hold(ax, "on"); grid(ax, "off"); box(ax, "on");
xlim(ax, xLimits); ylim(ax, yLimits);
khCurve = linspace(max(0.05, 0.02 .* xLimits(2)), xLimits(2), 600);
legendHandles = gobjects(1, numel(epsilonValues));
legendLabels = strings(1, numel(epsilonValues));
for ee = 1:numel(epsilonValues)
    style = cfg.curveStyles(min(ee, numel(cfg.curveStyles)));
    lineWidth = 1.1; if ee == 1, lineWidth = 1.5; end
    legendHandles(ee) = plot(ax, khCurve, predicted_ratio(khCurve, epsilonValues(ee)), style, ...
        "Color", cfg.curveColour, "LineWidth", lineWidth);
    legendLabels(ee) = sprintf("\\epsilon = %g%%", 100 .* epsilonValues(ee));
end
draw_xy_iqr(ax, kh, ratio, khErrLow, khErrHigh, ratioErrLow, ratioErrHigh, cfg.errorColour, yLimits);
scatter(ax, kh, ratio, cfg.markerArea, q, "filled", ...
    "MarkerEdgeColor", [0.15 0.15 0.15], "LineWidth", 0.5, "HandleVisibility", "off");
clim(ax, qLimits);
legendHandle = legend(ax, legendHandles, cellstr(legendLabels), ...
    "Location", "none", "FontSize", cfg.legendFontSize, "Box", "off", "Interpreter", "tex");
title(legendHandle, "Predicted", "FontSize", cfg.legendFontSize, "FontWeight", "normal");
xlabel(ax, "$kh=2\pi h/\lambda$", "Interpreter", "latex", "FontSize", cfg.labelFontSize);
ylabel(ax, "$h_{\mathrm{est}}/h$", "Interpreter", "latex", "FontSize", cfg.labelFontSize);
add_panel_label(ax, "(a)", cfg.panelLabelFontSize, [0.03 0.97], "top");
format_axes(ax, cfg);
end

function plot_sensitivity_panel(ax, khValues, referenceEpsilon, cfg)
epsilonCurve = linspace(0, 0.15, 400);
hold(ax, "on"); grid(ax, "off"); box(ax, "on");
xlim(ax, [0 15]); ylim(ax, [0 1.05]);
if referenceEpsilon > 0
    plot(ax, 100 .* [referenceEpsilon referenceEpsilon], [0 1.05], ":", ...
        "Color", cfg.errorColour, "LineWidth", 0.8, "HandleVisibility", "off");
    text(ax, 100 .* referenceEpsilon + 0.4, 1.02, sprintf("%g%%", 100 .* referenceEpsilon), ...
        "FontSize", cfg.annotationFontSize, "Color", [0.35 0.35 0.35], ...
        "HorizontalAlignment", "left", "VerticalAlignment", "top", "Interpreter", "tex");
end
for kk = 1:numel(khValues)
    style = cfg.curveStyles(min(kk, numel(cfg.curveStyles)));
    r = predicted_ratio(khValues(kk), epsilonCurve);
    plot(ax, 100 .* epsilonCurve, r, style, "Color", cfg.curveColour, ...
        "LineWidth", 1.1, "HandleVisibility", "off");
    labelEpsilon = 0.78 .* epsilonCurve(end);
    labelRatio = predicted_ratio(khValues(kk), labelEpsilon) + 0.04;
    text(ax, 100 .* labelEpsilon, labelRatio, sprintf("{\\itkh} = %g", khValues(kk)), ...
        "FontSize", cfg.annotationFontSize, "Color", cfg.curveColour, ...
        "HorizontalAlignment", "left", "VerticalAlignment", "bottom", "Interpreter", "tex");
end
xlabel(ax, "$\epsilon$ (\%)", "Interpreter", "latex", "FontSize", cfg.labelFontSize);
ylabel(ax, "$h_{\mathrm{est}}/h$", "Interpreter", "latex", "FontSize", cfg.labelFontSize);
add_panel_label(ax, "(b)", cfg.panelLabelFontSize, [0.97 0.97], "top", "right");
format_axes(ax, cfg);
end

function position_curve_legend(ax, legendHandle)
drawnow;
oldAxesUnits = ax.Units; oldLegendUnits = legendHandle.Units;
ax.Units = "normalized"; legendHandle.Units = "normalized";
axPosition = ax.Position; legendPosition = legendHandle.Position;
legendPosition(1) = axPosition(1) + axPosition(3) - legendPosition(3) - 0.04 .* axPosition(3);
legendPosition(2) = axPosition(2) + axPosition(4) - legendPosition(4) - 0.05 .* axPosition(4);
legendHandle.Position = legendPosition;
ax.Units = oldAxesUnits; legendHandle.Units = oldLegendUnits;
end

function draw_xy_iqr(ax, x, y, xLow, xHigh, yLow, yHigh, errorColour, limits)
if isempty(x), return; end
capHalfWidth = 0.012 .* diff(xlim(ax));
capHalfHeight = 0.012 .* diff(limits);
for ii = 1:numel(x)
    xEnds = [x(ii)-xLow(ii), x(ii)+xHigh(ii)];
    yEnds = [y(ii)-yLow(ii), y(ii)+yHigh(ii)];
    plot(ax, xEnds, [y(ii) y(ii)], "Color", errorColour, "LineWidth", 0.7, "HandleVisibility", "off");
    plot(ax, [xEnds(1) xEnds(1)], [y(ii)-capHalfHeight y(ii)+capHalfHeight], "Color", errorColour, "LineWidth", 0.7, "HandleVisibility", "off");
    plot(ax, [xEnds(2) xEnds(2)], [y(ii)-capHalfHeight y(ii)+capHalfHeight], "Color", errorColour, "LineWidth", 0.7, "HandleVisibility", "off");
    plot(ax, [x(ii) x(ii)], yEnds, "Color", errorColour, "LineWidth", 0.7, "HandleVisibility", "off");
    xCap = [x(ii)-capHalfWidth x(ii)+capHalfWidth];
    plot(ax, xCap, [yEnds(1) yEnds(1)], "Color", errorColour, "LineWidth", 0.7, "HandleVisibility", "off");
    plot(ax, xCap, [yEnds(2) yEnds(2)], "Color", errorColour, "LineWidth", 0.7, "HandleVisibility", "off");
end
end

function limits = padded_limits(lowerValues, upperValues)
lowerValues = lowerValues(isfinite(lowerValues)); upperValues = upperValues(isfinite(upperValues));
if isempty(lowerValues) || isempty(upperValues), limits = [0 1]; return; end
lower = min(lowerValues); upper = max(upperValues);
if upper <= lower, upper = lower + 1; end
padding = 0.07 .* (upper-lower);
limits = [lower-padding upper+padding];
end

function add_panel_label(ax, labelText, fontSize, position, verticalAlignment, horizontalAlignment)
if nargin < 6 || strlength(string(horizontalAlignment)) == 0, horizontalAlignment = "left"; end
text(ax, position(1), position(2), labelText, "Units", "normalized", ...
    "HorizontalAlignment", horizontalAlignment, "VerticalAlignment", verticalAlignment, ...
    "FontSize", fontSize, "FontWeight", "bold", "Interpreter", "none", "Clipping", "off");
end

function format_axes(ax, cfg)
axis(ax, "square");
set(ax, "FontSize", cfg.axesFontSize, "LineWidth", 0.8, "TickDir", "out", ...
    "TickLabelInterpreter", "tex", "Layer", "top");
end

function map = local_coolwarm(m)
if nargin < 1 || isempty(m), m = 256; end
anchors = [59 76 192; 84 112 222; 129 164 251; 180 205 251; ...
    221 221 221; 241 184 156; 229 112 88; 203 62 56; 180 4 38] ./ 255;
map = interp1(linspace(0, 1, size(anchors, 1)), anchors, linspace(0, 1, m), "linear");
map = max(0, min(1, map));
end
