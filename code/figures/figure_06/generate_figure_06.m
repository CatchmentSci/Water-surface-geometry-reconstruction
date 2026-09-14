function outputs = generate_figure_06(dataRoot, outputDir)
%GENERATE_FIGURE_06 Reproduce the six-panel accepted real-case WSG figure.
% Create a six-panel figure of accepted real-world WSE maps for:
%   R1, R3, R5, R7, R9 and R11.
%
% Case numbering matches calling_real_batch_analysis_latex_table.m:
%   1) sort the lookup by decreasing discharge;
%   2) retain Q < 124;
%   3) renumber the retained subset as R1-R13.
%
% Each panel shows the accepted WSE as a change in elevation relative to
% the flat initial WSE, following panel b of
% KLT_plot_refined_wave_solution_three_panel_batch_consistent.m:
%
%   delta WSE = accepted WSE - median flat-WSE level
%
% All panels use the same coolwarm colormap and the same colour limits,
% calculated from the finite minimum and maximum across all six maps. The
% original manually selected spatial focus and panel positions are specified
% explicitly below. Each checkpoint is transformed into the fixed R1
% streamwise/cross-stream frame so identical x/y coordinates represent the
% same physical position in every panel. Equal physical scaling is retained.
% Profile/transect lines are deliberately not drawn.
%
% The figure is exported at the approximately 5.5-inch text width of
% \documentclass[draft]{agujournal2019}. Its original 8.5:11 aspect ratio
% is retained, so all normalized panel geometry and map extents are unchanged.
%
% Outputs:
%   real_accepted_wse_R1_R3_R5_R7_R9_R11.png
%   real_accepted_wse_R1_R3_R5_R7_R9_R11.pdf
%
% Required companion files:
%   KLT_select_wave_checkpoint_refined.m
%   klt_analysis_case_lookup.tsv

%% ------------------------------------------------------------------------
% Inputs and repository paths
%% ------------------------------------------------------------------------

arguments
    dataRoot (1,1) string
    outputDir (1,1) string = string(fullfile(fileparts(mfilename('fullpath')), 'output'))
end

inputFolder = char(fullfile(dataRoot, 'videos', 'outputs'));
outputFolder = char(outputDir);

thisScriptDir = fileparts(mfilename('fullpath'));
lookupFile = char(fullfile(dataRoot, 'videos', 'inputs', ...
    'klt_analysis_case_lookup.tsv'));

checkpointPattern = '*_checkpoint.mat';
maxDischargeExclusive = 124;
requestedCaseNumbers = [1 3 5 7 9 11];
commonReferenceCaseNumber = 1;

outputDpi = 600;
axisTickSpacing_m = 5;

% Manually adjusted figure geometry, original local-coordinate focus windows
% and colour-bar position. The focus windows are transformed into the common
% R1-aligned frame below, so the same WSE regions remain displayed.
% Rows correspond to R1, R3, R5, R7, R9 and R11, respectively.
manualFigurePositionIn = [ ...
    0.5, 0.5, 5.5, 5.5 .* 11 ./ 8.5];

manualXLim = [ ...
     2.87878742561843, 24.5685182561474; ...
     5.38799004093502, 28.9601414776259; ...
     9.15823519021163, 26.1449166606539; ...
     7.13735647626544, 30.4153874166438; ...
     5.07880622754865, 27.7102251973610; ...
     8.06638105056105, 31.9418202168462];

manualYLim = [ ...
    -25.8601644924947,  -6.33888169706823; ...
    -24.3803437847717,  -3.34426383199207; ...
    -27.0385561067801, -11.8794311323655; ...
    -28.3596744486097,  -8.77567444860966; ...
    -31.6354069800956, -12.5954069800956; ...
    -31.8502007463063, -11.7635962637311];


manualCommonXLim = [2.87878742561851 24.5685182561477;3.07854465983637 23.6322589334301;5.78076808892576 21.3795929897661;2.42036059216062 23.9372047383662;2.38704234923035 23.9909214367994;3.50133365594289 22.8301977560394];
manualCommonYLim = [-25.8601644924946 -6.33888169706761;-25.5526903847382 -6.52799324324689;-25.1764044343067 -10.8591980446196;-24.7196537467289 -5.71995325534385;-24.3260714439246 -6.01438661973731;-23.8666076221467 -6.8251306680902];
%manualAxesPosition = [0.0912308341261861 0.7173190203072 0.352538331747628 0.2451809796928;0.571103044679204 0.7173190203072 0.342793910641592 0.2451809796928;0.0946520881172122 0.427574912672844 0.345695823765575 0.245180979692801;0.562835464419197 0.427574912672844 0.359329071161606 0.245180979692801;0.0803309771543711 0.137830805038488 0.374338045691257 0.2451809796928;0.562559202542203 0.137830805038488 0.359881594915593 0.2451809796928];


% The middle and bottom rows are raised by 0.018 and 0.036, respectively,
% to reserve sufficient space for the bottom x-axis labels and colourbar.
manualAxesPosition = [ ...
    0.0912308504848256, 0.717319020307200, 0.352538299030349, 0.245180979692800; ...
    0.5647273354745150, 0.717319020307200, 0.355545329050969, 0.245180979692800; ...
    0.0897273732987923, 0.427574912672844, 0.355545253402415, 0.245180979692801; ...
    0.5539287880762510, 0.427574912672844, 0.377142423847498, 0.245180979692801; ...
    0.0789287880762508, 0.137830805038488, 0.377142423847498, 0.245180979692800; ...
    0.5539288069883890, 0.137830805038488, 0.377142386023221, 0.245180979692800];

manualColorbarPosition = [ ...
    0.0789287880762509, 0.0500000000000000, ...
    0.8521424238474980, 0.0151515151515152];

% Original figure-normalized positions of the three manually added regions.
% These are converted below into data-linked polygons in the common frame.
% The R5 and R9 positions include the same vertical offsets applied to their
% axes above.
rectangleCaseNumbers = [5 3 9];
manualRectanglePosition = [ ...
    0.143156862745098, 0.559666666666670, ...
    0.122774509803921, 0.0407196969696971; ...
    0.730166666666667, 0.784090909090909, ...
    0.0921372549019608, 0.0388257575757576; ...
    0.220588235294118, 0.184674242424242, ...
    0.102941176470588, 0.0454545454545455];

rectangleEdgeColor = [0 0 0];
rectangleLineWidth = 0.5;
rectangleLineStyle = '-';

outputBase = fullfile(outputFolder, 'Figure6');
outputPngFile = [outputBase, '.png'];
outputPdfFile = [outputBase, '.pdf'];

% Display controls copied from panel b of the supplied template.
maskAcceptedWseForDisplay = true;
usePhysicalAxesForMapPlots = true;
plotColormap = local_coolwarm(256);

% Consistent typography and axes styling.
fontName = 'Arial';
axesFontSize = 9;
labelFontSize = 10;
panelLabelFontSize = 10;
axesLineWidth = 0.75;

%% ------------------------------------------------------------------------
% Selector configuration
%% ------------------------------------------------------------------------
if ~isempty(thisScriptDir)
    addpath(thisScriptDir, '-begin');
end
clear KLT_select_wave_checkpoint_refined
rehash;

selectorPath = which('KLT_select_wave_checkpoint_refined');
if isempty(selectorPath)
    error(['KLT_select_wave_checkpoint_refined.m was not found. Place it ', ...
           'beside this script or add its folder to the MATLAB path.']);
end
fprintf('\nUsing selector: %s\n', selectorPath);

% Keep these identical to calling_real_batch_analysis_latex_table.m.
selectorOpts = struct;
selectorOpts.stabilityAssessmentMode = ...
    'first_populated_10m_from_camera';
selectorOpts.firstSensedLength_m = 10;
selectorOpts.cameraRowOrder = 'ascending';
selectorOpts.displaySummary = false;
selectorOpts.usePmusicWavelength = true;
selectorOpts.pmusicModelOrder = 2;
selectorOpts.pmusicNfft = 2048;
selectorOpts.pmusicWindowLength = [];
selectorOpts.pmusicOverlap = [];
selectorOpts.requireForwardAmplitudePlateau = true;
selectorOpts.weakForwardWindow = 5;
selectorOpts.weakMaxForwardAmpGrowth = 0.05;
selectorOpts.weakMaxForwardAmpAbsGrowth_m = 0.003;
selectorOpts.strongForwardWindow = 5;
selectorOpts.strongMaxForwardAmpGrowth = 0.05;
selectorOpts.strongMaxForwardAmpAbsGrowth_m = 0.003;
selectorOpts.noPlateauFallbackMode = 'latest_supported';

%% ------------------------------------------------------------------------
% Validate paths and reproduce the real-case numbering
%% ------------------------------------------------------------------------
if ~isfolder(inputFolder)
    error('Input folder does not exist: %s', inputFolder);
end
if ~isfile(lookupFile)
    error('Lookup table does not exist: %s', lookupFile);
end
if ~isfolder(outputFolder)
    mkdir(outputFolder);
end

lookup = readtable(lookupFile, ...
    'FileType', 'text', ...
    'Delimiter', '\t', ...
    'CommentStyle', '#', ...
    'TextType', 'string', ...
    'VariableNamingRule', 'preserve');

requiredLookupColumns = ["filename_in", "sweep_value"];
missingColumns = setdiff(requiredLookupColumns, ...
    string(lookup.Properties.VariableNames));
if ~isempty(missingColumns)
    error('Lookup table is missing required column(s): %s', ...
        char(strjoin(missingColumns, ', ')));
end

lookup.filename_in = strtrim(string(lookup.filename_in));
if isnumeric(lookup.sweep_value)
    lookup.sweep_value = double(lookup.sweep_value);
else
    lookup.sweep_value = str2double(string(lookup.sweep_value));
end

blankRows = strlength(lookup.filename_in) == 0 & ...
    ~isfinite(lookup.sweep_value);
lookup(blankRows, :) = [];

if isempty(lookup)
    error('The lookup table contains no data rows.');
end
if any(strlength(lookup.filename_in) == 0)
    error('The lookup table contains a blank filename_in value.');
end
if any(~isfinite(lookup.sweep_value))
    error('Every lookup sweep_value must be a finite numeric value.');
end

lookup.lookupKey = strings(height(lookup), 1);
for ii = 1:height(lookup)
    lookup.lookupKey(ii) = ...
        normalise_checkpoint_key(lookup.filename_in(ii));
end

if numel(unique(lower(lookup.lookupKey))) ~= height(lookup)
    error('The lookup table contains duplicate filename_in values.');
end
if numel(unique(lookup.sweep_value)) ~= height(lookup)
    error(['The lookup table contains duplicate sweep_value values. ', ...
           'Unique discharges are required for an unambiguous R ranking.']);
end

lookup = sortrows(lookup, 'sweep_value', 'descend');
lookup = lookup(lookup.sweep_value < maxDischargeExclusive, :);
lookup.caseNumber = (1:height(lookup)).';
lookup.caseLabel = "R" + string(lookup.caseNumber);

if any(~ismember(requestedCaseNumbers, lookup.caseNumber))
    missingCaseNumbers = requestedCaseNumbers( ...
        ~ismember(requestedCaseNumbers, lookup.caseNumber));
    error('Requested case number(s) are unavailable after filtering: %s', ...
        char(strjoin(string(missingCaseNumbers), ', ')));
end
if requestedCaseNumbers(1) ~= commonReferenceCaseNumber
    error(['The common reference case must be the first entry in ', ...
           'requestedCaseNumbers so its coordinate frame is loaded first.']);
end

%% ------------------------------------------------------------------------
% Match lookup rows to checkpoint files
%% ------------------------------------------------------------------------
checkpointFiles = dir(fullfile(inputFolder, checkpointPattern));
checkpointFiles = checkpointFiles(~[checkpointFiles.isdir]);
if isempty(checkpointFiles)
    error('No files matched %s in: %s', checkpointPattern, inputFolder);
end

fileKeys = strings(numel(checkpointFiles), 1);
for ii = 1:numel(checkpointFiles)
    fileKeys(ii) = normalise_checkpoint_key(checkpointFiles(ii).name);
end

if numel(unique(lower(fileKeys))) ~= numel(fileKeys)
    error(['Two or more checkpoint files reduce to the same lookup key. ', ...
           'Resolve the duplicate filenames before plotting.']);
end

requestedRows = NaN(size(requestedCaseNumbers));
requestedFilePaths = strings(size(requestedCaseNumbers));
for ii = 1:numel(requestedCaseNumbers)
    rowIdx = find(lookup.caseNumber == requestedCaseNumbers(ii), 1, 'first');
    requestedRows(ii) = rowIdx;

    fileIdx = find(strcmpi(lookup.lookupKey(rowIdx), fileKeys), 1, 'first');
    if isempty(fileIdx)
        error('No checkpoint file was found for %s (%s).', ...
            char(lookup.caseLabel(rowIdx)), ...
            char(lookup.filename_in(rowIdx)));
    end

    requestedFilePaths(ii) = string(fullfile( ...
        checkpointFiles(fileIdx).folder, checkpointFiles(fileIdx).name));
end

%% ------------------------------------------------------------------------
% Load accepted maps before plotting so global limits can be calculated
%% ------------------------------------------------------------------------
nCases = numel(requestedCaseNumbers);
if ~isequal(size(manualXLim), [nCases, 2]) || ...
        ~isequal(size(manualYLim), [nCases, 2]) || ...
        ~isequal(size(manualAxesPosition), [nCases, 4])
    error(['Manual limits and positions must contain one row for each ', ...
           'requested case.']);
end
if size(manualRectanglePosition, 2) ~= 4 || ...
        size(manualRectanglePosition, 1) ~= numel(rectangleCaseNumbers)
    error('Each manual rectangle position must contain [x y width height].');
end
if any(~ismember(rectangleCaseNumbers, requestedCaseNumbers))
    error('Every rectangleCaseNumbers entry must be a requested case.');
end
caseData = repmat(empty_case_data(), nCases, 1);
allFiniteWse = [];
commonOrigin = [];
commonEs = [];
commonEn = [];

for ii = 1:nCases
    caseNumber = requestedCaseNumbers(ii);
    lookupRow = requestedRows(ii);
    caseLabel = lookup.caseLabel(lookupRow);
    filePath = char(requestedFilePaths(ii));

    fprintf('\n============================================================\n');
    fprintf('Preparing %s: %s\n', char(caseLabel), filePath);

    vars = who('-file', filePath);
    if ~ismember('checkpoint', vars)
        error('Checkpoint variable is missing from: %s', filePath);
    end

    S = load(filePath, 'checkpoint');
    [result, mapTable] = KLT_select_wave_checkpoint_refined( ...
        S.checkpoint, selectorOpts);

    [stopMap, selectedMapRow] = ...
        resolve_stop_map_for_plot(result, mapTable);

    if isfield(S.checkpoint, 'aa') && isfinite(S.checkpoint.aa)
        aa = round(S.checkpoint.aa);
    else
        aa = 1;
    end

    if ~isfield(S.checkpoint, 'wse_map') || ...
            size(S.checkpoint.wse_map, 1) < aa || ...
            size(S.checkpoint.wse_map, 2) < stopMap
        error('checkpoint.wse_map{%d,%d} is unavailable for %s.', ...
            aa, stopMap, char(caseLabel));
    end

    Waccepted = double(S.checkpoint.wse_map{aa, stopMap});
    [Ny, Nx] = size(Waccepted);

    requiredGeometryFields = {'origin','e_s','e_n','X_rot','Y_rot','xi','yi'};
    missingGeometryFields = requiredGeometryFields( ...
        ~isfield(S.checkpoint, requiredGeometryFields));
    if ~isempty(missingGeometryFields)
        error('Checkpoint geometry is incomplete for %s: %s', ...
            char(caseLabel), strjoin(missingGeometryFields, ', '));
    end
    if ~isequal(size(S.checkpoint.X_rot), size(Waccepted)) || ...
            ~isequal(size(S.checkpoint.Y_rot), size(Waccepted)) || ...
            numel(S.checkpoint.xi) ~= Nx || numel(S.checkpoint.yi) ~= Ny
        error('Checkpoint coordinate dimensions do not match the WSE map for %s.', ...
            char(caseLabel));
    end

    caseOrigin = double(S.checkpoint.origin(:).');
    caseEs = double(S.checkpoint.e_s(:).');
    caseEn = double(S.checkpoint.e_n(:).');
    validate_coordinate_basis(caseOrigin, caseEs, caseEn, char(caseLabel));

    if caseNumber == commonReferenceCaseNumber
        commonOrigin = caseOrigin;
        commonEs = caseEs;
        commonEn = caseEn;
        fprintf('  Common coordinate reference established from %s.\n', ...
            char(caseLabel));
    elseif isempty(commonOrigin)
        error('The common coordinate reference has not yet been loaded.');
    end

    if isfield(result, 'referenceMask') && ...
            isequal(size(result.referenceMask), size(Waccepted))
        referenceMask = logical(result.referenceMask);
    else
        referenceMask = isfinite(Waccepted);
    end

    flatWseLevel = get_flat_wse_level( ...
        S.checkpoint, aa, referenceMask, Waccepted);
    WacceptedRelative = Waccepted - flatWseLevel;
    WacceptedDisplay = WacceptedRelative;

    if maskAcceptedWseForDisplay
        WacceptedDisplay(~referenceMask) = NaN;
    end

    % Preserve the existing colour-limit calculation before cropping the map
    % to its original panel focus window.
    finiteVals = WacceptedDisplay(isfinite(WacceptedDisplay));
    if isempty(finiteVals)
        error('Accepted WSE display map contains no finite values for %s.', ...
            char(caseLabel));
    end

    if ~usePhysicalAxesForMapPlots
        error(['Common spatial registration requires ', ...
               'usePhysicalAxesForMapPlots = true.']);
    end

    [xMap, yMap] = world_grid_to_common_frame( ...
        S.checkpoint.X_rot, S.checkpoint.Y_rot, ...
        commonOrigin, commonEs, commonEn);

    % Hide values outside the original local-coordinate viewing rectangle.
    % This retains precisely the WSE region that was visible before the
    % common-frame transformation.
    [localXGrid, localYGrid] = meshgrid( ...
        double(S.checkpoint.xi(:).'), double(S.checkpoint.yi(:)));
    originalFocusMask = ...
        localXGrid >= manualXLim(ii, 1) & ...
        localXGrid <= manualXLim(ii, 2) & ...
        localYGrid >= manualYLim(ii, 1) & ...
        localYGrid <= manualYLim(ii, 2);
    WacceptedPlot = WacceptedDisplay;
    WacceptedPlot(~originalFocusMask) = NaN;

    [commonXLim, commonYLim] = local_limits_to_common_frame( ...
        manualXLim(ii, :), manualYLim(ii, :), ...
        caseOrigin, caseEs, caseEn, ...
        commonOrigin, commonEs, commonEn);

    % Apply the manually refined viewport in the common R1 coordinate frame.
    commonXLim = manualCommonXLim(ii, :);
    commonYLim = manualCommonYLim(ii, :);

    caseData(ii).caseNumber = caseNumber;
    caseData(ii).caseLabel = caseLabel;
    caseData(ii).Q = lookup.sweep_value(lookupRow);
    caseData(ii).filePath = string(filePath);
    caseData(ii).selectedMapIndex = stopMap;
    caseData(ii).selectedMapRow = selectedMapRow;
    caseData(ii).stopMethod = string(result.stopMethod);
    caseData(ii).flatWseLevel_m = flatWseLevel;
    caseData(ii).WacceptedRelative_m = WacceptedPlot;
    caseData(ii).xMap = xMap;
    caseData(ii).yMap = yMap;
    caseData(ii).xLabel = "$x$ (m)";
    caseData(ii).yLabel = "$y$ (m)";
    caseData(ii).xLimits = commonXLim;
    caseData(ii).yLimits = commonYLim;
    caseData(ii).localOrigin = caseOrigin;
    caseData(ii).localEs = caseEs;
    caseData(ii).localEn = caseEn;

    allFiniteWse = [allFiniteWse; finiteVals(:)]; %#ok<AGROW>

    fprintf('  Accepted map index: %d (map-table row %d)\n', ...
        stopMap, selectedMapRow);
    fprintf('  Stop method: %s\n', char(caseData(ii).stopMethod));
    fprintf('  Flat WSE level: %.6g m\n', flatWseLevel);
    fprintf('  Relative WSE range: %.6g to %.6g m\n', ...
        min(finiteVals), max(finiteVals));
    fprintf('  Common-frame panel x limits: %.6g to %.6g\n', ...
        caseData(ii).xLimits(1), caseData(ii).xLimits(2));
    fprintf('  Common-frame panel y limits: %.6g to %.6g\n', ...
        caseData(ii).yLimits(1), caseData(ii).yLimits(2));
end

globalColorLimits = [min(allFiniteWse), max(allFiniteWse)];
if ~all(isfinite(globalColorLimits))
    error('The combined WSE colour limits are not finite.');
end
if globalColorLimits(1) == globalColorLimits(2)
    pad = max(abs(globalColorLimits(1)) .* 0.05, 1e-6);
    globalColorLimits = globalColorLimits + [-pad, pad];
end

fprintf('\nShared colour limits: %.6g to %.6g m\n', ...
    globalColorLimits(1), globalColorLimits(2));
fprintf(['Each panel retains its original local-coordinate focus window, ', ...
         'transformed into the common R%d frame.\n\n'], ...
    commonReferenceCaseNumber);

% Adjust panel widths to the transformed limit aspect ratios while retaining
% the existing panel centres, heights and row positions.
plotAxesPosition = manualAxesPosition;
figureWidthIn = manualFigurePositionIn(3);
figureHeightIn = manualFigurePositionIn(4);
for ii = 1:nCases
    xSpan = diff(caseData(ii).xLimits);
    ySpan = diff(caseData(ii).yLimits);
    newWidth = plotAxesPosition(ii, 4) .* ...
        (figureHeightIn ./ figureWidthIn) .* (xSpan ./ ySpan);
    oldCentre = manualAxesPosition(ii, 1) + manualAxesPosition(ii, 3) ./ 2;
    plotAxesPosition(ii, 1) = oldCentre - newWidth ./ 2;
    plotAxesPosition(ii, 3) = newWidth;
end
if any(plotAxesPosition(:, 1) < 0) || ...
        any(sum(plotAxesPosition(:, [1 3]), 2) > 1)
    error('A transformed panel position extends beyond the figure canvas.');
end

% Convert the existing figure-normalized region boxes into data-linked
% polygons in the common frame. They therefore remain attached to the same
% WSE features when a case is rotated into the reference orientation.
nRectangles = size(manualRectanglePosition, 1);
rectangleCommonX = cell(nRectangles, 1);
rectangleCommonY = cell(nRectangles, 1);
for rr = 1:nRectangles
    caseIndex = find( ...
        requestedCaseNumbers == rectangleCaseNumbers(rr), 1, 'first');
    [localRectX, localRectY] = annotation_box_to_local_polygon( ...
        manualRectanglePosition(rr, :), ...
        manualAxesPosition(caseIndex, :), ...
        manualXLim(caseIndex, :), manualYLim(caseIndex, :));
    [rectangleCommonX{rr}, rectangleCommonY{rr}] = ...
        local_points_to_common_frame( ...
            localRectX, localRectY, ...
            caseData(caseIndex).localOrigin, ...
            caseData(caseIndex).localEs, ...
            caseData(caseIndex).localEn, ...
            commonOrigin, commonEs, commonEn);
end

%% ------------------------------------------------------------------------
% Create the 3 x 2 AGU-page-aspect figure
%% ------------------------------------------------------------------------
fig = figure('Color', 'w', ...
    'Units', 'inches', ...
    'Position', manualFigurePositionIn, ...
    'PaperUnits', 'inches', ...
    'PaperPosition', [0 0 manualFigurePositionIn(3:4)], ...
    'PaperSize', manualFigurePositionIn(3:4), ...
    'InvertHardcopy', 'off');

axesHandles = gobjects(nCases, 1);
for ii = 1:nCases
    ax = axes('Parent', fig, ...
        'Units', 'normalized', ...
        'Position', plotAxesPosition(ii, :), ...
        'PositionConstraint', 'innerposition');
    axesHandles(ii) = ax;

    hMap = surface(ax, ...
        caseData(ii).xMap, ...
        caseData(ii).yMap, ...
        zeros(size(caseData(ii).WacceptedRelative_m)), ...
        caseData(ii).WacceptedRelative_m, ...
        'EdgeColor', 'none', ...
        'FaceColor', 'texturemap');
    view(ax, 2);
    hold(ax, 'on');
    colormap(ax, plotColormap);
    set(hMap, ...
        'AlphaData', isfinite(caseData(ii).WacceptedRelative_m), ...
        'FaceAlpha', 'texturemap', ...
        'AlphaDataMapping', 'none');

    rectangleRows = find(rectangleCaseNumbers == caseData(ii).caseNumber);
    for rr = rectangleRows(:).'
        plot(ax, rectangleCommonX{rr}, rectangleCommonY{rr}, ...
            'Color', rectangleEdgeColor, ...
            'LineWidth', rectangleLineWidth, ...
            'LineStyle', rectangleLineStyle);
    end

    set(ax, ...
        'YDir', 'normal', ...
        'Color', 'w', ...
        'FontName', fontName, ...
        'FontSize', axesFontSize, ...
        'LineWidth', axesLineWidth, ...
        'Box', 'on', ...
        'TickDir', 'out', ...
        'Layer', 'top');

    clim(ax, globalColorLimits);
    daspect(ax, [1 1 1]);

    % Set the limits after the data aspect ratio so they remain exact.
    xlim(ax, caseData(ii).xLimits);
    ylim(ax, caseData(ii).yLimits);
    set(ax, ...
        'XLimMode', 'manual', ...
        'YLimMode', 'manual', ...
        'XTick', ticks_at_fixed_spacing( ...
            caseData(ii).xLimits, axisTickSpacing_m), ...
        'YTick', ticks_at_fixed_spacing( ...
            caseData(ii).yLimits, axisTickSpacing_m), ...
        'XTickMode', 'manual', ...
        'YTickMode', 'manual');
    xtickformat(ax, '%.0f');
    ytickformat(ax, '%.0f');

    add_panel_label(ax, sprintf('(%c)', 'a' + ii - 1), ...
        fontName, panelLabelFontSize);

    rowNumber = ceil(ii ./ 2);
    columnNumber = mod(ii - 1, 2) + 1;

    % Retain tick values on every panel, but avoid repeated axis titles.
    if rowNumber == 3
        xlabel(ax, char(caseData(ii).xLabel), ...
            'Interpreter', 'latex', ...
            'FontName', fontName, ...
            'FontSize', labelFontSize);
    end
    if columnNumber == 1
        ylabel(ax, char(caseData(ii).yLabel), ...
            'Interpreter', 'latex', ...
            'FontName', fontName, ...
            'FontSize', labelFontSize);
    end
end

% One horizontal colour ramp beneath all six panels.
cb = colorbar(axesHandles(end), 'southoutside');
cb.Units = 'normalized';
cb.Position = manualColorbarPosition;
cb.FontName = fontName;
cb.FontSize = axesFontSize;
cb.Label.String = '$\Delta z_{\mathrm{WSE}}$ (m)';
cb.Label.FontName = fontName;
cb.Label.FontSize = labelFontSize;
cb.Label.Interpreter = 'latex';

% Creating a colour bar can alter its associated axes. Reapply every
% manually measured position after the colour bar has been created.
drawnow;
for ii = 1:nCases
    axesHandles(ii).Units = 'normalized';
    axesHandles(ii).PositionConstraint = 'innerposition';
    axesHandles(ii).Position = plotAxesPosition(ii, :);
end
cb.Units = 'normalized';
cb.Position = manualColorbarPosition;

%% ------------------------------------------------------------------------
% Save outputs
%% ------------------------------------------------------------------------
exportgraphics(fig, outputPngFile, 'Resolution', outputDpi);
exportgraphics(fig, outputPdfFile, 'ContentType', 'vector');

fprintf('Six-panel figure written to:\n');
fprintf('  %s\n', outputPngFile);
fprintf('  %s\n', outputPdfFile);

outputs = struct( ...
    'pngFile', string(outputPngFile), ...
    'pdfFile', string(outputPdfFile), ...
    'lookupFile', string(lookupFile), ...
    'caseNumbers', requestedCaseNumbers, ...
    'checkpointFiles', requestedFilePaths, ...
    'discharges_m3s', [caseData.Q], ...
    'selectedMapIndices', [caseData.selectedMapIndex]);
close(fig);
end

%% ========================================================================
% Local helper functions
%% ========================================================================
function D = empty_case_data()
    D = struct;
    D.caseNumber = NaN;
    D.caseLabel = "";
    D.Q = NaN;
    D.filePath = "";
    D.selectedMapIndex = NaN;
    D.selectedMapRow = NaN;
    D.stopMethod = "";
    D.flatWseLevel_m = NaN;
    D.WacceptedRelative_m = [];
    D.xMap = [];
    D.yMap = [];
    D.xLabel = "";
    D.yLabel = "";
    D.xLimits = [NaN NaN];
    D.yLimits = [NaN NaN];
    D.localOrigin = [NaN NaN];
    D.localEs = [NaN NaN];
    D.localEn = [NaN NaN];
end

function key = normalise_checkpoint_key(fileName)
    [~, key, ~] = fileparts(char(strtrim(string(fileName))));
    key = regexprep(key, '_checkpoint$', '', 'ignorecase');
    key = regexprep(key, '_ckpt$', '', 'ignorecase');
    key = string(strtrim(key));
end

function [stopMap, selectedMapRow] = ...
        resolve_stop_map_for_plot(result, mapTable)
    selectedMapRow = NaN;

    if isfield(result, 'selectedMapRow') && ...
            isfinite(result.selectedMapRow)
        selectedMapRow = round(result.selectedMapRow);
    end

    if ~isfield(result, 'selectedMapIndex') || ...
            ~isfinite(result.selectedMapIndex)
        error('result.selectedMapIndex is missing or non-finite.');
    end
    stopMap = round(result.selectedMapIndex);

    if ~isfinite(selectedMapRow) && istable(mapTable) && ...
            ismember('mapIndex', mapTable.Properties.VariableNames)
        selectedMapRow = find( ...
            round(mapTable.mapIndex) == stopMap, 1, 'first');
    end

    if ~isfinite(selectedMapRow)
        error('The selected map could not be resolved to a map-table row.');
    end
end

function flatWseLevel = get_flat_wse_level( ...
        checkpoint, aa, referenceMask, Waccepted)
    flatWseLevel = NaN;

    try
        if isfield(checkpoint, 'wse_map') && ...
                size(checkpoint.wse_map, 1) >= aa && ...
                size(checkpoint.wse_map, 2) >= 1
            W0 = checkpoint.wse_map{aa, 1};
            if isnumeric(W0) && ~isempty(W0) && ...
                    isequal(size(W0), size(Waccepted))
                vals = double(W0(referenceMask & isfinite(W0)));
                if isempty(vals)
                    vals = double(W0(isfinite(W0)));
                end
                if ~isempty(vals)
                    flatWseLevel = median(vals, 'omitnan');
                end
            end
        end
    catch
    end

    if ~isfinite(flatWseLevel)
        vals = Waccepted(referenceMask & isfinite(Waccepted));
        if isempty(vals)
            vals = Waccepted(isfinite(Waccepted));
        end
        if isempty(vals)
            flatWseLevel = 0;
        else
            flatWseLevel = median(vals, 'omitnan');
        end
    end
end

function validate_coordinate_basis(origin, eS, eN, caseLabel)
    if numel(origin) ~= 2 || numel(eS) ~= 2 || numel(eN) ~= 2 || ...
            any(~isfinite([origin eS eN]))
        error('Invalid two-dimensional coordinate basis for %s.', caseLabel);
    end
    if abs(norm(eS) - 1) > 1e-6 || abs(norm(eN) - 1) > 1e-6 || ...
            abs(dot(eS, eN)) > 1e-6
        error('Coordinate basis is not orthonormal for %s.', caseLabel);
    end
end

function [xCommon, yCommon] = world_grid_to_common_frame( ...
        Xworld, Yworld, commonOrigin, commonEs, commonEn)
    dX = double(Xworld) - commonOrigin(1);
    dY = double(Yworld) - commonOrigin(2);
    xCommon = dX .* commonEs(1) + dY .* commonEs(2);
    yCommon = dX .* commonEn(1) + dY .* commonEn(2);
end

function [xCommon, yCommon] = local_points_to_common_frame( ...
        xLocal, yLocal, caseOrigin, caseEs, caseEn, ...
        commonOrigin, commonEs, commonEn)
    Xworld = caseOrigin(1) + ...
        xLocal .* caseEs(1) + yLocal .* caseEn(1);
    Yworld = caseOrigin(2) + ...
        xLocal .* caseEs(2) + yLocal .* caseEn(2);
    [xCommon, yCommon] = world_grid_to_common_frame( ...
        Xworld, Yworld, commonOrigin, commonEs, commonEn);
end

function [commonXLim, commonYLim] = local_limits_to_common_frame( ...
        localXLim, localYLim, caseOrigin, caseEs, caseEn, ...
        commonOrigin, commonEs, commonEn)
    xCorners = [localXLim(1), localXLim(2), ...
                localXLim(2), localXLim(1)];
    yCorners = [localYLim(1), localYLim(1), ...
                localYLim(2), localYLim(2)];
    [xCommon, yCommon] = local_points_to_common_frame( ...
        xCorners, yCorners, caseOrigin, caseEs, caseEn, ...
        commonOrigin, commonEs, commonEn);
    commonXLim = [min(xCommon), max(xCommon)];
    commonYLim = [min(yCommon), max(yCommon)];
end

function [xLocal, yLocal] = annotation_box_to_local_polygon( ...
        boxPosition, axesPosition, localXLim, localYLim)
    figureX = boxPosition(1) + ...
        [0, boxPosition(3), boxPosition(3), 0, 0];
    figureY = boxPosition(2) + ...
        [0, 0, boxPosition(4), boxPosition(4), 0];
    axesX = (figureX - axesPosition(1)) ./ axesPosition(3);
    axesY = (figureY - axesPosition(2)) ./ axesPosition(4);
    xLocal = localXLim(1) + axesX .* diff(localXLim);
    yLocal = localYLim(1) + axesY .* diff(localYLim);
end

function ticks = ticks_at_fixed_spacing(limits, spacing)
    firstTick = ceil(limits(1) ./ spacing) .* spacing;
    lastTick = floor(limits(2) ./ spacing) .* spacing;
    ticks = firstTick:spacing:lastTick;
end

function add_panel_label(ax, txt, fontName, fontSize)
    % Match the notation and inside-top-left placement of the companion plot.
    text(ax, 0.02, 0.97, txt, ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', ...
        'FontName', fontName, ...
        'FontSize', fontSize, ...
        'FontWeight', 'bold', ...
        'Interpreter', 'none', ...
        'Clipping', 'on');
end

function map = local_coolwarm(m)
    if nargin < 1 || isempty(m)
        m = 256;
    end
    if m <= 0
        map = zeros(0, 3);
        return
    end

    % Exact compact coolwarm anchors used by the supplied panel-b template.
    anchors = [ ...
        59   76  192; ...
        84  112  222; ...
        129 164  251; ...
        180 205  251; ...
        221 221  221; ...
        241 184  156; ...
        229 112   88; ...
        203  62   56; ...
        180   4   38] ./ 255;

    x = linspace(0, 1, size(anchors, 1));
    xi = linspace(0, 1, m);
    map = interp1(x, anchors, xi, 'linear');
    map = max(0, min(1, map));
end
