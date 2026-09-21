function outputs = reproduce_hydraulic_statistics(archiveRoot, outputFolder)
%REPRODUCE_HYDRAULIC_STATISTICS Rebuild the hydraulic input to Table 3.
%   OUTPUTS = REPRODUCE_HYDRAULIC_STATISTICS(ARCHIVEROOT, OUTPUTFOLDER)
%   reads the deposited observation inventory, stage/discharge workbook,
%   cross-section survey and long-term discharge values.
%
% For each immediate subfolder whose name contains a timestamp formatted as
% yyyymmdd_HHMM (for example, "Video 2 - 20180315_0818"), this script:
%   1. finds the adjacent workbook timestamps bracketing each video timestamp;
%   2. uses an exact workbook record where available, otherwise linearly
%      interpolates Stage, Local survey stage and Discharge;
%   3. calculates cross-sectional flow geometry from the surveyed profile;
%   4. calculates the flow exceedance percentage from the long-term flow record;
%   5. calculates mean velocity, Froude number and Reynolds numbers; and
%   6. writes the results to XLSX and CSV files.
%
% Cross-section area, water-surface top width and wetted perimeter are
% calculated by linearly interpolating the waterline crossing within each
% surveyed cross-section segment.
%
% Reynolds-number note:
%   The original file described Re = U * R_h / nu but calculated using
%   hydraulic depth. This version reports both:
%       ReynoldsNumber_Rh  = U * R_h / nu
%       ReynoldsNumber_4Rh = U * (4 * R_h) / nu
%   The second is a commonly used open-channel convention and is used for
%   the Reynolds-regime label below.

arguments
    archiveRoot (1, 1) string
    outputFolder (1, 1) string = fullfile(fileparts(mfilename('fullpath')), "output")
end

inputFolder = fullfile(archiveRoot, "videos", "inputs");
hydraulicFolder = fullfile(inputFolder, "hydraulics");
flowDataFile = fullfile(hydraulicFolder, "Austins_SG_Q.xlsx");
crossSectionFile = fullfile(inputFolder, "cross_section.csv");
longTermFlowDataFile = fullfile(hydraulicFolder, "long_term_discharge_values.csv");
observationInventoryFile = fullfile(hydraulicFolder, "video_observation_inventory.csv");
if ~isfolder(outputFolder)
    mkdir(outputFolder);
end
outputCsvFile = fullfile(outputFolder, "Dart_video_hydraulic_statistics.csv");

% Settings.
waterTemp_C = 10;              % degrees Celsius
maxMatchOffset_minutes = 30;   % flag if either interpolation bracket is farther away

% Optional quality filter for the long-term flow record. The default empty
% value uses every finite, non-negative discharge in the attached CSV. To
% exclude provisional values, for example, use ["Good"; "Estimated"].
longTermAcceptedQuality = strings(0, 1);

% Constants.
g = 9.81;                      % gravitational acceleration, m/s^2
nu = water_kinematic_viscosity(waterTemp_C);

%% ------------------------------------------------------------------------
% CHECK INPUT FILES AND FOLDERS
% -------------------------------------------------------------------------

if exist(flowDataFile, 'file') ~= 2
    error('Stage/discharge workbook not found: %s', flowDataFile);
end
if exist(crossSectionFile, 'file') ~= 2
    error('Cross-section CSV not found. Checked: %s', crossSectionFile);
end
if exist(longTermFlowDataFile, 'file') ~= 2
    error('Long-term flow CSV not found. Checked: %s', longTermFlowDataFile);
end
if exist(observationInventoryFile, 'file') ~= 2
    error('Video observation inventory not found: %s', observationInventoryFile);
end

%% ------------------------------------------------------------------------
% READ STAGE AND DISCHARGE DATA
% -------------------------------------------------------------------------

flowTable = readtable(flowDataFile, 'VariableNamingRule', 'preserve');
if width(flowTable) < 4
    error('The workbook must contain at least four columns: Datetime, Stage, Local survey stage and Discharge.');
end

flowDatetime = parse_datetime_column(flowTable{:, 1});
stage_m = parse_numeric_column(flowTable{:, 2});
localSurveyStage_m = parse_numeric_column(flowTable{:, 3});
discharge_m3s = parse_numeric_column(flowTable{:, 4});

validFlowRows = ~isnat(flowDatetime) & isfinite(stage_m) & ...
    isfinite(localSurveyStage_m) & isfinite(discharge_m3s);

if ~all(validFlowRows)
    warning('Discarding %d workbook rows with missing or invalid values.', sum(~validFlowRows));
end

flowDatetime = flowDatetime(validFlowRows);
stage_m = stage_m(validFlowRows);
localSurveyStage_m = localSurveyStage_m(validFlowRows);
discharge_m3s = discharge_m3s(validFlowRows);

[flowDatetime, flowSortOrder] = sort(flowDatetime);
stage_m = stage_m(flowSortOrder);
localSurveyStage_m = localSurveyStage_m(flowSortOrder);
discharge_m3s = discharge_m3s(flowSortOrder);

if isempty(flowDatetime)
    error('No usable timestamped stage/discharge records were found in the workbook.');
end

%% ------------------------------------------------------------------------
% READ LONG-TERM FLOW RECORD FOR EXCEEDANCE CALCULATIONS
% -------------------------------------------------------------------------

% For a matched discharge Q, flow exceedance is calculated as:
%   100 * number of long-term observations >= Q / number of observations
% This is the empirical exceedance probability used in a flow-duration curve.
longTermFlowTable = readtable(longTermFlowDataFile, 'VariableNamingRule', 'preserve');
longTermVariableNames = string(longTermFlowTable.Properties.VariableNames);

longTermValueColumn = find(strcmpi(longTermVariableNames, "value"), 1, 'first');
if isempty(longTermValueColumn)
    error('The long-term flow CSV must contain a discharge column named value.');
end

longTermDischarge_m3s = parse_numeric_column(longTermFlowTable{:, longTermValueColumn});
validLongTermRows = isfinite(longTermDischarge_m3s) & longTermDischarge_m3s >= 0;

longTermQualityColumn = find(strcmpi(longTermVariableNames, "quality"), 1, 'first');
if ~isempty(longTermAcceptedQuality) && ~isempty(longTermQualityColumn)
    longTermQuality = strtrim(string(longTermFlowTable{:, longTermQualityColumn}));
    validLongTermRows = validLongTermRows & ...
        ismember(lower(longTermQuality), lower(longTermAcceptedQuality));
end

longTermDischarge_m3s = longTermDischarge_m3s(validLongTermRows);
if isempty(longTermDischarge_m3s)
    error('No usable discharge values were found in the long-term flow CSV.');
end

fprintf('Using %d long-term flow observations to calculate exceedance percentages.\n', ...
    numel(longTermDischarge_m3s));
if ~isempty(longTermAcceptedQuality) && ~isempty(longTermQualityColumn)
    fprintf('Retained long-term quality flags: %s\n', ...
        char(strjoin(longTermAcceptedQuality, ', ')));
end

%% ------------------------------------------------------------------------
% READ CROSS-SECTION SURVEY
% -------------------------------------------------------------------------

crossSectionTable = readtable(crossSectionFile, 'VariableNamingRule', 'preserve');
if width(crossSectionTable) < 3
    error('The cross-section CSV must contain x, y and elevation columns.');
end

crossSection_x = parse_numeric_column(crossSectionTable{:, 1});
crossSection_y = parse_numeric_column(crossSectionTable{:, 2});
crossSectionElevation_m = parse_numeric_column(crossSectionTable{:, 3});

validCrossSectionRows = isfinite(crossSection_x) & isfinite(crossSection_y) & ...
    isfinite(crossSectionElevation_m);

crossSection_x = crossSection_x(validCrossSectionRows);
crossSection_y = crossSection_y(validCrossSectionRows);
crossSectionElevation_m = crossSectionElevation_m(validCrossSectionRows);

if numel(crossSection_x) < 2
    error('The cross-section CSV must contain at least two valid surveyed points.');
end

% Cumulative horizontal distance along the surveyed cross-section line.
profileSegmentLength_m = hypot(diff(crossSection_x), diff(crossSection_y));
if any(profileSegmentLength_m <= 0)
    warning('The cross-section contains repeated consecutive x/y points. Zero-length segments will be ignored.');
end
profileStation_m = [0; cumsum(profileSegmentLength_m)];

%% ------------------------------------------------------------------------
% READ VIDEO OBSERVATION INVENTORY
% -------------------------------------------------------------------------
observationInventory = readtable(observationInventoryFile, ...
    'VariableNamingRule', 'preserve', 'TextType', 'string');
requiredInventoryFields = ["FolderName", "VideoDatetime"];
missingInventoryFields = setdiff(requiredInventoryFields, ...
    string(observationInventory.Properties.VariableNames));
if ~isempty(missingInventoryFields)
    error('Observation inventory is missing field(s): %s', ...
        strjoin(missingInventoryFields, ', '));
end
folderNames = string(observationInventory.FolderName);
folderDatetime = parse_datetime_column(observationInventory.VideoDatetime);
if any(strlength(folderNames) == 0 | isnat(folderDatetime))
    error('Observation inventory contains a blank folder name or invalid datetime.');
end

[folderDatetime, folderSortOrder] = sort(folderDatetime);
folderNames = folderNames(folderSortOrder);

%% ------------------------------------------------------------------------
% MATCH EACH VIDEO FOLDER AND CALCULATE HYDRAULIC STATISTICS
% -------------------------------------------------------------------------

nVideos = numel(folderDatetime);

% Matching and interpolation audit fields.
% MatchedFlowDatetime and MatchOffset_minutes retain the nearest supporting
% workbook record for convenient quality checking. The hydraulic values are
% calculated from exact values or bracketed linear interpolation.
matchedFlowDatetime = NaT(nVideos, 1);
beforeFlowDatetime = NaT(nVideos, 1);
afterFlowDatetime = NaT(nVideos, 1);
timeOffset_minutes = nan(nVideos, 1);
beforeOffset_minutes = nan(nVideos, 1);
afterOffset_minutes = nan(nVideos, 1);
bracketSpan_minutes = nan(nVideos, 1);
interpolationFraction = nan(nVideos, 1);
interpolationUsed = false(nVideos, 1);
matchedWithinTolerance = false(nVideos, 1);
dataStatus = strings(nVideos, 1);

matchedStage_m = nan(nVideos, 1);
matchedLocalSurveyStage_m = nan(nVideos, 1);
matchedDischarge_m3s = nan(nVideos, 1);
flowExceedance_percent = nan(nVideos, 1);

meanVelocity_mps = nan(nVideos, 1);
flowArea_m2 = nan(nVideos, 1);
topWidth_m = nan(nVideos, 1);
averageDepth_m = nan(nVideos, 1);
wettedPerimeter_m = nan(nVideos, 1);
hydraulicRadius_m = nan(nVideos, 1);
FroudeNumber = nan(nVideos, 1);
ReynoldsNumber_Rh = nan(nVideos, 1);
ReynoldsNumber_4Rh = nan(nVideos, 1);

for iVideo = 1:nVideos
    videoTime = folderDatetime(iVideo);

    % Use the workbook values directly when the timestamps are synchronous.
    exactRow = find(flowDatetime == videoTime, 1, 'first');

    if ~isempty(exactRow)
        beforeRow = exactRow;
        afterRow = exactRow;

        matchedFlowDatetime(iVideo) = flowDatetime(exactRow);
        beforeFlowDatetime(iVideo) = flowDatetime(exactRow);
        afterFlowDatetime(iVideo) = flowDatetime(exactRow);
        timeOffset_minutes(iVideo) = 0;
        beforeOffset_minutes(iVideo) = 0;
        afterOffset_minutes(iVideo) = 0;
        bracketSpan_minutes(iVideo) = 0;
        interpolationFraction(iVideo) = 0;
        interpolationUsed(iVideo) = false;
        matchedWithinTolerance(iVideo) = true;
        dataStatus(iVideo) = "Exact workbook timestamp";

        matchedStage_m(iVideo) = stage_m(exactRow);
        matchedLocalSurveyStage_m(iVideo) = localSurveyStage_m(exactRow);
        matchedDischarge_m3s(iVideo) = discharge_m3s(exactRow);

    else
        % For a non-synchronous timestamp, find the immediately preceding
        % and following workbook records. Do not extrapolate beyond the
        % available workbook time range.
        beforeRow = find(flowDatetime < videoTime, 1, 'last');
        afterRow = find(flowDatetime > videoTime, 1, 'first');

        if isempty(beforeRow) || isempty(afterRow)
            dataStatus(iVideo) = "Outside workbook timestamp range";
            continue;
        end

        beforeFlowDatetime(iVideo) = flowDatetime(beforeRow);
        afterFlowDatetime(iVideo) = flowDatetime(afterRow);
        beforeOffset_minutes(iVideo) = minutes(videoTime - flowDatetime(beforeRow));
        afterOffset_minutes(iVideo) = minutes(flowDatetime(afterRow) - videoTime);
        bracketSpan_minutes(iVideo) = minutes(flowDatetime(afterRow) - flowDatetime(beforeRow));
        interpolationFraction(iVideo) = beforeOffset_minutes(iVideo) / ...
            bracketSpan_minutes(iVideo);

        % Retain the closest supporting record for the existing match audit
        % fields, while calculating values from both bracketing records.
        [timeOffset_minutes(iVideo), nearestBracketIndex] = min( ...
            [beforeOffset_minutes(iVideo), afterOffset_minutes(iVideo)]);
        bracketRows = [beforeRow, afterRow];
        matchedFlowDatetime(iVideo) = flowDatetime(bracketRows(nearestBracketIndex));

        matchedWithinTolerance(iVideo) = max( ...
            beforeOffset_minutes(iVideo), afterOffset_minutes(iVideo)) <= ...
            maxMatchOffset_minutes;
        interpolationUsed(iVideo) = true;
        dataStatus(iVideo) = "Linearly interpolated";

        fraction = interpolationFraction(iVideo);
        matchedStage_m(iVideo) = stage_m(beforeRow) + ...
            fraction * (stage_m(afterRow) - stage_m(beforeRow));
        matchedLocalSurveyStage_m(iVideo) = localSurveyStage_m(beforeRow) + ...
            fraction * (localSurveyStage_m(afterRow) - localSurveyStage_m(beforeRow));
        matchedDischarge_m3s(iVideo) = discharge_m3s(beforeRow) + ...
            fraction * (discharge_m3s(afterRow) - discharge_m3s(beforeRow));
    end

    % All downstream statistics use the exact or interpolated values above.
    flowExceedance_percent(iVideo) = calculate_flow_exceedance_percent( ...
        longTermDischarge_m3s, matchedDischarge_m3s(iVideo));

    [flowArea_m2(iVideo), topWidth_m(iVideo), wettedPerimeter_m(iVideo)] = ...
        calculate_cross_section_geometry(profileStation_m, ...
        crossSectionElevation_m, matchedLocalSurveyStage_m(iVideo));

    if flowArea_m2(iVideo) > 0
        meanVelocity_mps(iVideo) = matchedDischarge_m3s(iVideo) / flowArea_m2(iVideo);
    end

    if topWidth_m(iVideo) > 0
        averageDepth_m(iVideo) = flowArea_m2(iVideo) / topWidth_m(iVideo);
    end

    if wettedPerimeter_m(iVideo) > 0
        hydraulicRadius_m(iVideo) = flowArea_m2(iVideo) / wettedPerimeter_m(iVideo);
    end

    if isfinite(meanVelocity_mps(iVideo)) && averageDepth_m(iVideo) > 0
        FroudeNumber(iVideo) = meanVelocity_mps(iVideo) / sqrt(g * averageDepth_m(iVideo));
    end

    if isfinite(meanVelocity_mps(iVideo)) && hydraulicRadius_m(iVideo) > 0
        ReynoldsNumber_Rh(iVideo) = meanVelocity_mps(iVideo) * hydraulicRadius_m(iVideo) / nu;
        ReynoldsNumber_4Rh(iVideo) = 4 * ReynoldsNumber_Rh(iVideo);
    end
end

%% ------------------------------------------------------------------------
% FLOW CLASSIFICATION
% -------------------------------------------------------------------------

FroudeRegime = classify_froude(FroudeNumber);
ReynoldsRegime_4Rh = classify_reynolds(ReynoldsNumber_4Rh);

%% ------------------------------------------------------------------------
% OUTPUT TABLE
% -------------------------------------------------------------------------

resultsTable = table( ...
    folderNames, ...
    folderDatetime, ...
    dataStatus, ...
    interpolationUsed, ...
    beforeFlowDatetime, ...
    afterFlowDatetime, ...
    interpolationFraction, ...
    bracketSpan_minutes, ...
    beforeOffset_minutes, ...
    afterOffset_minutes, ...
    matchedFlowDatetime, ...
    timeOffset_minutes, ...
    matchedWithinTolerance, ...
    matchedStage_m, ...
    matchedLocalSurveyStage_m, ...
    matchedDischarge_m3s, ...
    flowExceedance_percent, ...
    meanVelocity_mps, ...
    flowArea_m2, ...
    topWidth_m, ...
    averageDepth_m, ...
    wettedPerimeter_m, ...
    hydraulicRadius_m, ...
    FroudeNumber, ...
    ReynoldsNumber_Rh, ...
    ReynoldsNumber_4Rh, ...
    FroudeRegime, ...
    ReynoldsRegime_4Rh, ...
    'VariableNames', { ...
    'FolderName', ...
    'VideoDatetime', ...
    'DataStatus', ...
    'InterpolationUsed', ...
    'BeforeFlowDatetime', ...
    'AfterFlowDatetime', ...
    'InterpolationFraction', ...
    'BracketSpan_minutes', ...
    'BeforeOffset_minutes', ...
    'AfterOffset_minutes', ...
    'NearestFlowDatetime', ...
    'NearestRecordOffset_minutes', ...
    'InterpolationWithinTolerance', ...
    'Stage_m', ...
    'LocalSurveyStage_m', ...
    'Discharge_m3s', ...
    'FlowExceedance_percent', ...
    'MeanVelocity_mps', ...
    'FlowArea_m2', ...
    'TopWidth_m', ...
    'AverageDepth_m', ...
    'WettedPerimeter_m', ...
    'HydraulicRadius_m', ...
    'FroudeNumber', ...
    'ReynoldsNumber_Rh', ...
    'ReynoldsNumber_4Rh', ...
    'FroudeRegime', ...
    'ReynoldsRegime_4Rh'});

writetable(resultsTable, outputCsvFile);

fprintf('\nProcessed %d timestamped video folders.\n', nVideos);
fprintf('Wrote CSV results:  %s\n', outputCsvFile);

outsideRange = dataStatus == "Outside workbook timestamp range";
if any(outsideRange)
    warning('%d video folder timestamp(s) fall outside the workbook time range and were not extrapolated.', ...
        sum(outsideRange));
end

outsideTolerance = ~matchedWithinTolerance & ~outsideRange;
if any(outsideTolerance)
    warning('%d interpolated video folder(s) have at least one bracketing workbook timestamp more than %.1f minutes away. Check the interpolation audit columns.', ...
        sum(outsideTolerance), maxMatchOffset_minutes);
end

outputs = struct('csvFile', outputCsvFile, 'table', resultsTable, ...
    'flowDataFile', flowDataFile, ...
    'longTermFlowDataFile', longTermFlowDataFile, ...
    'crossSectionFile', crossSectionFile, ...
    'observationInventoryFile', observationInventoryFile);
end

%% ------------------------------------------------------------------------
% LOCAL FUNCTIONS
% -------------------------------------------------------------------------

function exceedancePercent = calculate_flow_exceedance_percent(longTermDischarge_m3s, discharge_m3s)
% Calculate the empirical percentage of long-term flow observations that
% are greater than or equal to a specified discharge.

    if ~isfinite(discharge_m3s) || isempty(longTermDischarge_m3s)
        exceedancePercent = NaN;
        return;
    end

    exceedancePercent = 100 * sum(longTermDischarge_m3s >= discharge_m3s) / ...
        numel(longTermDischarge_m3s);
end

function numericColumn = parse_numeric_column(rawColumn)
% Convert a table column to a numeric column vector.

    if isnumeric(rawColumn)
        numericColumn = double(rawColumn);
    else
        numericColumn = str2double(string(rawColumn));
    end
    numericColumn = numericColumn(:);
end

function datetimeColumn = parse_datetime_column(rawColumn)
% Convert common workbook timestamp representations to datetime values.

    if isdatetime(rawColumn)
        datetimeColumn = rawColumn;
        datetimeColumn = datetimeColumn(:);
        return;
    end

    if isnumeric(rawColumn)
        datetimeColumn = datetime(rawColumn, 'ConvertFrom', 'excel');
        datetimeColumn = datetimeColumn(:);
        return;
    end

    timestampText = string(rawColumn);
    timestampText = timestampText(:);
    datetimeColumn = NaT(size(timestampText));

    candidateFormats = { ...
        'yyyy-MM-dd''T''HH:mm:ss', ...
        'yyyy-MM-dd HH:mm:ss', ...
        'dd/MM/yyyy HH:mm:ss', ...
        'yyyy/MM/dd HH:mm:ss'};

    for iFormat = 1:numel(candidateFormats)
        unresolved = isnat(datetimeColumn) & ~ismissing(timestampText);
        if ~any(unresolved)
            break;
        end

        try
            datetimeColumn(unresolved) = datetime(timestampText(unresolved), ...
                'InputFormat', candidateFormats{iFormat});
        catch
            % Try the next explicit format.
        end
    end

    unresolved = isnat(datetimeColumn) & ~ismissing(timestampText);
    if any(unresolved)
        try
            datetimeColumn(unresolved) = datetime(timestampText(unresolved));
        catch
            % Invalid values remain NaT and will be discarded by the caller.
        end
    end
end

function [flowArea_m2, topWidth_m, wettedPerimeter_m] = ...
    calculate_cross_section_geometry(profileStation_m, bedElevation_m, waterLevel_m)
% Calculate wetted flow geometry for a piecewise-linear cross-section.
%
% flowArea_m2       = integrated submerged cross-sectional area
% topWidth_m        = summed water-surface width across submerged segments
% wettedPerimeter_m = submerged bed length

    flowArea_m2 = 0;
    topWidth_m = 0;
    wettedPerimeter_m = 0;

    for iSegment = 1:(numel(profileStation_m) - 1)
        segmentWidth_m = profileStation_m(iSegment + 1) - profileStation_m(iSegment);
        segmentRise_m = bedElevation_m(iSegment + 1) - bedElevation_m(iSegment);

        if segmentWidth_m <= 0
            continue;
        end

        depth1_m = waterLevel_m - bedElevation_m(iSegment);
        depth2_m = waterLevel_m - bedElevation_m(iSegment + 1);

        if depth1_m <= 0 && depth2_m <= 0
            % Entire segment is dry or exactly at the waterline.
            continue;

        elseif depth1_m >= 0 && depth2_m >= 0
            % Entire segment is submerged.
            flowArea_m2 = flowArea_m2 + 0.5 * (depth1_m + depth2_m) * segmentWidth_m;
            topWidth_m = topWidth_m + segmentWidth_m;
            wettedPerimeter_m = wettedPerimeter_m + hypot(segmentWidth_m, segmentRise_m);

        else
            % The waterline crosses this segment. Interpolate its location.
            crossingFraction = depth1_m / (depth1_m - depth2_m);

            if depth1_m > 0
                submergedWidth_m = segmentWidth_m * crossingFraction;
                submergedRise_m = segmentRise_m * crossingFraction;
                flowArea_m2 = flowArea_m2 + 0.5 * depth1_m * submergedWidth_m;
            else
                submergedWidth_m = segmentWidth_m * (1 - crossingFraction);
                submergedRise_m = segmentRise_m * (1 - crossingFraction);
                flowArea_m2 = flowArea_m2 + 0.5 * depth2_m * submergedWidth_m;
            end

            topWidth_m = topWidth_m + submergedWidth_m;
            wettedPerimeter_m = wettedPerimeter_m + hypot(submergedWidth_m, submergedRise_m);
        end
    end
end

function regime = classify_froude(froudeNumber)
% Classify flow using a +/-0.05 band around critical flow.

    regime = strings(size(froudeNumber));
    regime(:) = "Unavailable";
    regime(froudeNumber < 0.95) = "Subcritical";
    regime(abs(froudeNumber - 1) <= 0.05) = "Near critical";
    regime(froudeNumber > 1.05) = "Supercritical";
end

function regime = classify_reynolds(reynoldsNumber)
% Classify flow using the reported ReynoldsNumber_4Rh values.

    regime = strings(size(reynoldsNumber));
    regime(:) = "Unavailable";
    regime(reynoldsNumber < 500) = "Laminar";
    regime(reynoldsNumber >= 500 & reynoldsNumber < 2000) = "Transitional";
    regime(reynoldsNumber >= 2000) = "Turbulent";
end

function nu = water_kinematic_viscosity(tempC)
% Return approximate freshwater kinematic viscosity in m^2/s.

    T = tempC + 273.15;  % Kelvin

    % Dynamic viscosity of water, Pa s.
    mu = 2.414e-5 * 10^(247.8 / (T - 140));

    % Approximate freshwater density, kg/m^3.
    rho = 1000 * (1 - ((tempC + 288.9414) / ...
        (508929.2 * (tempC + 68.12963))) * (tempC - 3.9863)^2);

    % Kinematic viscosity, m^2/s.
    nu = mu / rho;
end
