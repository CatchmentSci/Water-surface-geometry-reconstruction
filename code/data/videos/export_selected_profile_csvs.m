function outputs = export_selected_profile_csvs(archiveRoot, auxiliaryFolder, outputFolder)
%EXPORT_SELECTED_PROFILE_CSVS Rebuild the 13 deposited profile CSV files.
arguments
    archiveRoot (1,1) string
    auxiliaryFolder (1,1) string
    outputFolder (1,1) string
end
% Export one CSV per real case containing the exact selected-map profile
% metrics used by KLT_select_wave_checkpoint_refined to create the batch
% summary consumed by calling_real_batch_analysis_latex_table.m.
%
% Each output row is one profile line from profileTable at the selector's
% accepted/selected map. The script also applies the selector's aggregation
% rule to amplitude98_2_m and pmusicWavelength_m:
%   1) keep finite values;
%   2) when at least minProfiles values exist, reject values farther than
%      3.5 times the scaled MAD from the initial median;
%   3) calculate the retained median, Q25, Q75 and retained count.
%
% Before writing each CSV, these calculations are checked against the
% selector's result fields. A mismatch stops the export rather than writing
% data that are not comparable with the batch summary.
%
% Cross-section/depth geometry and velocity-cell counts are joined by
% profileIndex from the existing first-10-m transect CSV for the same case.
% rowCoord and transectN_m must also agree within validationTolerance, so an
% incorrectly paired auxiliary file cannot be written silently.

%% ------------------------------------------------------------------------
% User settings
%% ------------------------------------------------------------------------
archiveRoot = string(java.io.File(archiveRoot).getCanonicalPath());
inputFolder = fullfile(archiveRoot, "videos", "checkpoints");
lookupFile = fullfile(archiveRoot, "videos", "inputs", ...
    "klt_analysis_case_lookup.tsv");

checkpointPattern = '*_checkpoint.mat';
casePrefix = 'R';
maxDischargeExclusive = 124;
aggregationMinProfiles = 8; % Selector default used by the batch summary.
validationTolerance = 1e-12;
writeOutputs = true;

auxiliaryTransectCsvSuffix = ...
    '_first10m_transect_depth_from_observed_velocity_with_deep_water_velocity_pmusic_depths.csv';

% Set WSE_PROFILE_EXPORT_DRY_RUN=1 in the environment to run every
% calculation and validation without creating or replacing CSV files.
if strcmpi(strtrim(getenv('WSE_PROFILE_EXPORT_DRY_RUN')), '1')
    writeOutputs = false;
end

%% ------------------------------------------------------------------------
% Resolve companion files and selector
%% ------------------------------------------------------------------------
thisScriptDir = fileparts(mfilename('fullpath'));

selectorCandidates = {fullfile(thisScriptDir, ...
    'KLT_select_wave_checkpoint_refined.m')};

selectorFile = '';
for ii = 1:numel(selectorCandidates)
    if isfile(selectorCandidates{ii})
        selectorFile = selectorCandidates{ii};
        break
    end
end
if isempty(selectorFile)
    error('KLT_select_wave_checkpoint_refined.m was not found in any configured location.');
end

addpath(fileparts(selectorFile), '-begin');
clear KLT_select_wave_checkpoint_refined
rehash;
selectorFileResolved = which('KLT_select_wave_checkpoint_refined');
if isempty(selectorFileResolved)
    error('KLT_select_wave_checkpoint_refined.m could not be placed on the MATLAB path.');
end
selectorFileUsed = "code/data/videos/KLT_select_wave_checkpoint_refined.m";

fprintf('\nSelector: %s\n', selectorFileUsed);
fprintf('Lookup  : %s\n', lookupFile);
fprintf('Inputs  : %s\n', inputFolder);
fprintf('Outputs : %s\n', outputFolder);
fprintf('Mode    : %s\n\n', ternary_text(writeOutputs, 'write CSV files', 'dry run; do not write CSV files'));

%% ------------------------------------------------------------------------
% Selector configuration: matched to calling_real_batch_analysis_latex_table
%% ------------------------------------------------------------------------
selectorOpts = struct;
selectorOpts.stabilityAssessmentMode = 'first_populated_10m_from_camera';
selectorOpts.firstSensedLength_m = 10;
selectorOpts.cameraRowOrder = 'ascending';
selectorOpts.displaySummary = false;
selectorOpts.usePmusicWavelength = true;
selectorOpts.pmusicModelOrder = 2;
selectorOpts.pmusicNfft = 2048;
selectorOpts.pmusicWindowLength = [];
selectorOpts.pmusicOverlap = [];
selectorOpts.minProfiles = aggregationMinProfiles;
selectorOpts.requireForwardAmplitudePlateau = true;
selectorOpts.weakForwardWindow = 5;
selectorOpts.weakMaxForwardAmpGrowth = 0.05;
selectorOpts.weakMaxForwardAmpAbsGrowth_m = 0.003;
selectorOpts.strongForwardWindow = 5;
selectorOpts.strongMaxForwardAmpGrowth = 0.05;
selectorOpts.strongMaxForwardAmpAbsGrowth_m = 0.003;
selectorOpts.noPlateauFallbackMode = 'latest_supported';

%% ------------------------------------------------------------------------
% Read and filter the case lookup exactly as in the batch summary
%% ------------------------------------------------------------------------
if ~isfolder(inputFolder)
    error('Input folder does not exist: %s', inputFolder);
end
if ~isfile(lookupFile)
    error('Lookup file does not exist: %s', lookupFile);
end

lookup = readtable(lookupFile, ...
    'FileType', 'text', ...
    'Delimiter', '\t', ...
    'CommentStyle', '#', ...
    'TextType', 'string', ...
    'VariableNamingRule', 'preserve');

requiredLookupColumns = ["filename_in", "sweep_value"];
missingLookupColumns = setdiff(requiredLookupColumns, ...
    string(lookup.Properties.VariableNames));
if ~isempty(missingLookupColumns)
    error('Lookup is missing required column(s): %s', ...
        char(strjoin(missingLookupColumns, ', ')));
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

if any(strlength(lookup.filename_in) == 0) || ...
        any(~isfinite(lookup.sweep_value))
    error('Every lookup row must have a filename_in and finite sweep_value.');
end

lookup.lookupKey = strings(height(lookup), 1);
for ii = 1:height(lookup)
    lookup.lookupKey(ii) = normalise_checkpoint_key(lookup.filename_in(ii));
end

if numel(unique(lower(lookup.lookupKey))) ~= height(lookup)
    error('The lookup contains duplicate filename keys.');
end
if numel(unique(lookup.sweep_value)) ~= height(lookup)
    error('The lookup contains duplicate sweep values.');
end

lookup = sortrows(lookup, 'sweep_value', 'descend');
lookup = lookup(lookup.sweep_value < maxDischargeExclusive, :);
lookup.caseNumber = (1:height(lookup)).';
lookup.caseLabel = string(casePrefix) + string(lookup.caseNumber);

checkpointFiles = dir(fullfile(inputFolder, checkpointPattern));
checkpointFiles = checkpointFiles(~[checkpointFiles.isdir]);
fileKeys = strings(numel(checkpointFiles), 1);
for ii = 1:numel(checkpointFiles)
    fileKeys(ii) = normalise_checkpoint_key(checkpointFiles(ii).name);
end

if numel(unique(lower(fileKeys))) ~= numel(fileKeys)
    error('Two or more checkpoint files reduce to the same lookup key.');
end

if writeOutputs && ~isfolder(outputFolder)
    mkdir(outputFolder);
end

%% ------------------------------------------------------------------------
% Export selected-map profile data for every retained case
%% ------------------------------------------------------------------------
nWritten = 0;
nValidated = 0;
nMissing = 0;

for ii = 1:height(lookup)
    matchIdx = find(strcmpi(lookup.lookupKey(ii), fileKeys), 1, 'first');
    if isempty(matchIdx)
        warning('No checkpoint matched %s (%s).', ...
            lookup.caseLabel(ii), lookup.filename_in(ii));
        nMissing = nMissing + 1;
        continue
    end

    checkpointFileName = string(checkpointFiles(matchIdx).name);
    checkpointFilePath = fullfile( ...
        checkpointFiles(matchIdx).folder, checkpointFiles(matchIdx).name);

    fprintf('Processing %s: %s, Q = %.10g\n', ...
        lookup.caseLabel(ii), checkpointFileName, lookup.sweep_value(ii));

    vars = who('-file', checkpointFilePath);
    if ~ismember('checkpoint', vars)
        warning('Skipping %s: MAT file contains no checkpoint variable.', ...
            checkpointFileName);
        continue
    end

    loaded = load(checkpointFilePath, 'checkpoint');
    checkpoint = loaded.checkpoint;
    [result, mapTable, profileTable] = ...
        KLT_select_wave_checkpoint_refined(checkpoint, selectorOpts);

    [selectedMapIndex, selectedMapRow] = ...
        resolve_selected_map(result, mapTable);
    selectedMask = round(profileTable.mapIndex) == selectedMapIndex;
    selectedProfiles = profileTable(selectedMask, :);
    if isempty(selectedProfiles)
        warning('Skipping %s: selected map %d has no profileTable rows.', ...
            checkpointFileName, selectedMapIndex);
        continue
    end
    selectedProfiles = sortrows(selectedProfiles, 'profileIndex');

    requiredProfileColumns = ["profileIndex", "amplitude98_2_m", ...
        "autocorrWavelength_m", "autocorrPeakR", ...
        "pmusicWavelength_m", "pmusicPeakFrequency_cpm", ...
        "pmusicPeakPower", "pmusicStatus"];
    missingProfileColumns = setdiff(requiredProfileColumns, ...
        string(selectedProfiles.Properties.VariableNames));
    if ~isempty(missingProfileColumns)
        error('profileTable is missing column(s): %s', ...
            char(strjoin(missingProfileColumns, ', ')));
    end

    amplitudeAgg = aggregate_like_selector( ...
        selectedProfiles.amplitude98_2_m, aggregationMinProfiles);
    pmusicAgg = aggregate_like_selector( ...
        selectedProfiles.pmusicWavelength_m, aggregationMinProfiles);

    validate_aggregate_against_result( ...
        amplitudeAgg, result.amplitude_m, result.amplitudeIQR_m, ...
        result.nAmplitudeProfiles, validationTolerance, ...
        checkpointFileName + " amplitude");
    validate_aggregate_against_result( ...
        pmusicAgg, result.pmusicWavelength_m, result.pmusicWavelengthIQR_m, ...
        result.nPmusicWavelengthProfiles, validationTolerance, ...
        checkpointFileName + " PMUSIC wavelength");

    nProfiles = height(selectedProfiles);
    profileIndex = round(selectedProfiles.profileIndex);
    rowCoord = nan(nProfiles, 1);
    transectN_m = nan(nProfiles, 1);

    if isfield(result, 'rowCoord') && ~isempty(result.rowCoord)
        resultRowCoord = result.rowCoord(:);
        validProfileIndex = profileIndex >= 1 & ...
            profileIndex <= numel(resultRowCoord);
        rowCoord(validProfileIndex) = ...
            resultRowCoord(profileIndex(validProfileIndex));

        if isfield(checkpoint, 'yi') && ~isempty(checkpoint.yi)
            yi = double(checkpoint.yi(:));
            transectN_m(validProfileIndex) = interp1( ...
                (1:numel(yi)).', yi, rowCoord(validProfileIndex), ...
                'linear', NaN);
        end
    end

    auxiliaryTransect = load_auxiliary_transect_columns( ...
        auxiliaryFolder, lookup.filename_in(ii), auxiliaryTransectCsvSuffix, ...
        profileIndex, rowCoord, transectN_m, validationTolerance);

    profilePmusicVelocity_mps = nan(nProfiles, 1);
    validPmusic = isfinite(selectedProfiles.pmusicWavelength_m) & ...
        selectedProfiles.pmusicWavelength_m > 0;
    profilePmusicVelocity_mps(validPmusic) = sqrt( ...
        9.81 .* selectedProfiles.pmusicWavelength_m(validPmusic) ./ (2*pi));

    summaryUsDeepWater_mps = NaN;
    if isfinite(result.pmusicWavelength_m) && result.pmusicWavelength_m > 0
        summaryUsDeepWater_mps = sqrt( ...
            9.81 .* result.pmusicWavelength_m ./ (2*pi));
    end

    selectionReason = describe_selected_map_reason( ...
        result, selectedMapIndex, selectedMapRow);
    fallbackSelection = is_fallback_selection( ...
        result.stopMethod, selectionReason);

    amplitudeRetentionStatus = retention_status( ...
        selectedProfiles.amplitude98_2_m, amplitudeAgg.keep);
    pmusicRetentionStatus = retention_status( ...
        selectedProfiles.pmusicWavelength_m, pmusicAgg.keep);

    T = table( ...
        repmat(lookup.caseNumber(ii), nProfiles, 1), ...
        repmat(lookup.caseLabel(ii), nProfiles, 1), ...
        repmat(lookup.filename_in(ii), nProfiles, 1), ...
        repmat(lookup.sweep_value(ii), nProfiles, 1), ...
        repmat(checkpointFileName, nProfiles, 1), ...
        repmat(string(selectorFileUsed), nProfiles, 1), ...
        repmat(selectedMapIndex, nProfiles, 1), ...
        repmat(selectedMapRow, nProfiles, 1), ...
        repmat(string(result.stopMethod), nProfiles, 1), ...
        repmat(selectionReason, nProfiles, 1), ...
        repmat(fallbackSelection, nProfiles, 1), ...
        profileIndex, rowCoord, auxiliaryTransect.rowCell, transectN_m, ...
        auxiliaryTransect.nearestRowIndex, ...
        auxiliaryTransect.crossSectionX_m, ...
        auxiliaryTransect.crossSectionY_m, ...
        auxiliaryTransect.crossSectionBedElevation_m, ...
        auxiliaryTransect.crossSectionInitialMedianWSE_m, ...
        auxiliaryTransect.crossSectionDepth_m, ...
        auxiliaryTransect.crossSectionDepthStatus, ...
        auxiliaryTransect.nIntersectedCells, ...
        auxiliaryTransect.nFiniteVelocityCells, ...
        selectedProfiles.amplitude98_2_m, ...
        amplitudeAgg.keep, amplitudeRetentionStatus, ...
        selectedProfiles.autocorrWavelength_m, ...
        selectedProfiles.autocorrPeakR, ...
        selectedProfiles.pmusicWavelength_m, ...
        selectedProfiles.pmusicPeakFrequency_cpm, ...
        selectedProfiles.pmusicPeakPower, ...
        string(selectedProfiles.pmusicStatus), ...
        pmusicAgg.keep, pmusicRetentionStatus, ...
        profilePmusicVelocity_mps, ...
        repmat(result.amplitude_m, nProfiles, 1), ...
        repmat(result.amplitudeIQR_m(1), nProfiles, 1), ...
        repmat(result.amplitudeIQR_m(2), nProfiles, 1), ...
        repmat(result.nAmplitudeProfiles, nProfiles, 1), ...
        repmat(result.pmusicWavelength_m, nProfiles, 1), ...
        repmat(result.pmusicWavelengthIQR_m(1), nProfiles, 1), ...
        repmat(result.pmusicWavelengthIQR_m(2), nProfiles, 1), ...
        repmat(result.nPmusicWavelengthProfiles, nProfiles, 1), ...
        repmat(summaryUsDeepWater_mps, nProfiles, 1), ...
        repmat("matches selector result", nProfiles, 1), ...
        'VariableNames', { ...
        'caseNumber','caseLabel','filenameLookup','sweepValue_Q', ...
        'checkpointFileName','selectorFileUsed', ...
        'selectedMapIndex','selectedMapRow','stopMethod', ...
        'selectionReason','fallbackSelection', ...
        'profileIndex','rowCoord','rowCell','transectN_m', ...
        'nearestRowIndex','crossSectionX_m','crossSectionY_m', ...
        'crossSectionBedElevation_m','crossSectionInitialMedianWSE_m', ...
        'crossSectionDepth_m','crossSectionDepthStatus', ...
        'nIntersectedCells','nFiniteVelocityCells', ...
        'amplitude98_2_m','amplitudeRetainedForSummary', ...
        'amplitudeRetentionStatus', ...
        'autocorrWavelength_m','autocorrPeakR', ...
        'pmusicWavelength_m','pmusicPeakFrequency_cpm', ...
        'pmusicPeakPower','pmusicStatus', ...
        'pmusicRetainedForSummary','pmusicRetentionStatus', ...
        'velocityFromProfilePmusic_mps', ...
        'summaryAmplitudeMedian_m','summaryAmplitudeQ25_m', ...
        'summaryAmplitudeQ75_m','summaryNAmplitudeProfiles', ...
        'summaryPmusicWavelengthMedian_m', ...
        'summaryPmusicWavelengthQ25_m', ...
        'summaryPmusicWavelengthQ75_m', ...
        'summaryNPmusicProfiles','summaryUsDeepWater_mps', ...
        'summaryValidationStatus'});

    outputFileName = lookup.filename_in(ii) + ...
        "_selected_map_profiles_for_real_batch_summary.csv";
    outputFilePath = fullfile(outputFolder, outputFileName);

    nValidated = nValidated + 1;
    if writeOutputs
        writetable(T, outputFilePath);
        nWritten = nWritten + 1;
        fprintf('  wrote %d profiles: %s\n', nProfiles, outputFilePath);
    else
        fprintf('  validated %d profiles; would write: %s\n', ...
            nProfiles, outputFilePath);
    end
end

fprintf('\n============================================================\n');
fprintf('Cases selected by Q filter : %d\n', height(lookup));
fprintf('Cases validated            : %d\n', nValidated);
fprintf('Missing checkpoints        : %d\n', nMissing);
fprintf('CSV files written          : %d\n', nWritten);
if ~writeOutputs
    fprintf('Dry run complete: no CSV files were written.\n');
end
fprintf('============================================================\n');

if nMissing ~= 0 || nValidated ~= height(lookup)
    error('Selected-profile export did not complete all %d cases.', height(lookup));
end
outputs = struct('outputFolder', char(outputFolder), ...
    'nCases', nValidated, 'nFilesWritten', nWritten);
end

%% ========================================================================
% Local functions
%% ========================================================================
function key = normalise_checkpoint_key(fileName)
    [~, key, ~] = fileparts(char(strtrim(string(fileName))));
    key = regexprep(key, '_checkpoint$', '', 'ignorecase');
    key = regexprep(key, '_ckpt$', '', 'ignorecase');
    key = string(strtrim(key));
end

function [selectedMapIndex, selectedMapRow] = ...
        resolve_selected_map(result, mapTable)
    selectedMapIndex = NaN;
    selectedMapRow = NaN;
    if isfield(result, 'selectedMapIndex') && ...
            isfinite(result.selectedMapIndex)
        selectedMapIndex = round(result.selectedMapIndex);
    end
    if isfield(result, 'selectedMapRow') && ...
            isfinite(result.selectedMapRow)
        selectedMapRow = round(result.selectedMapRow);
    end
    if ~isfinite(selectedMapRow) && isfinite(selectedMapIndex) && ...
            istable(mapTable) && ...
            ismember('mapIndex', mapTable.Properties.VariableNames)
        selectedMapRow = find( ...
            round(mapTable.mapIndex) == selectedMapIndex, 1, 'first');
    end
    if ~isfinite(selectedMapIndex) || ~isfinite(selectedMapRow)
        error('The selector returned no resolvable selected map.');
    end
end

function selectionReason = describe_selected_map_reason( ...
        result, selectedMapIndex, selectedMapRow)
    selectionReason = "";
    if isfield(result, 'stopInfo') && isstruct(result.stopInfo) && ...
            isfield(result.stopInfo, 'selectionReason')
        selectionReason = strtrim(string(result.stopInfo.selectionReason));
    end
    if strlength(selectionReason) == 0 && isfield(result, 'stopMethod')
        selectionReason = replace(string(result.stopMethod), "_", " ");
    end
    if strlength(selectionReason) == 0
        selectionReason = "selected map reason unavailable";
    end
    selectionReason = selectionReason + " (map " + ...
        string(selectedMapIndex) + ", row " + string(selectedMapRow) + ")";
end

function tf = is_fallback_selection(stopMethod, selectionReason)
    stopMethod = lower(string(stopMethod));
    selectionReason = lower(string(selectionReason));
    tf = any(contains(stopMethod, "latest_supported")) || ...
         any(contains(stopMethod, "fallback")) || ...
         any(contains(selectionReason, "latest quality-supported")) || ...
         any(contains(selectionReason, "latest finite map")) || ...
         any(contains(selectionReason, "fallback"));
end

function A = aggregate_like_selector(values, minProfiles)
    values = double(values(:));
    finiteMask = isfinite(values);
    finiteValues = values(finiteMask);

    A.keep = false(size(values));
    A.median = NaN;
    A.q25 = NaN;
    A.q75 = NaN;
    A.n = 0;
    A.initialMedian = NaN;
    A.scaledMAD = NaN;
    A.filterApplied = false;

    if isempty(finiteValues)
        return
    end

    A.initialMedian = median(finiteValues, 'omitnan');
    A.scaledMAD = 1.4826 .* median( ...
        abs(finiteValues - A.initialMedian), 'omitnan');
    keepFinite = true(size(finiteValues));

    if numel(finiteValues) >= minProfiles && ...
            isfinite(A.scaledMAD) && A.scaledMAD > 0
        keepFinite = abs(finiteValues - A.initialMedian) <= ...
            3.5 .* A.scaledMAD;
        A.filterApplied = true;
    end

    retainedValues = finiteValues(keepFinite);
    finiteIndices = find(finiteMask);
    A.keep(finiteIndices(keepFinite)) = true;
    A.n = numel(retainedValues);

    if ~isempty(retainedValues)
        A.median = median(retainedValues, 'omitnan');
        A.q25 = prctile(retainedValues, 25);
        A.q75 = prctile(retainedValues, 75);
    end
end

function status = retention_status(values, keep)
    values = double(values(:));
    keep = logical(keep(:));
    status = repmat("retained for summary", size(values));
    status(~isfinite(values)) = "not retained: non-finite value";
    status(isfinite(values) & ~keep) = ...
        "not retained: outside 3.5 scaled-MAD threshold";
end

function validate_aggregate_against_result( ...
        aggregate, expectedMedian, expectedIQR, expectedCount, ...
        tolerance, label)
    validate_numeric_scalar(aggregate.median, expectedMedian, ...
        tolerance, label + " median");
    validate_numeric_scalar(aggregate.q25, expectedIQR(1), ...
        tolerance, label + " Q25");
    validate_numeric_scalar(aggregate.q75, expectedIQR(2), ...
        tolerance, label + " Q75");
    if aggregate.n ~= expectedCount
        error('%s count mismatch: reconstructed %d, selector result %d.', ...
            label, aggregate.n, expectedCount);
    end
end

function validate_numeric_scalar(actual, expected, tolerance, label)
    if isnan(actual) && isnan(expected)
        return
    end
    scale = max([1, abs(actual), abs(expected)]);
    if ~isfinite(actual) || ~isfinite(expected) || ...
            abs(actual - expected) > tolerance .* scale
        error('%s mismatch: reconstructed %.17g, selector result %.17g.', ...
            label, actual, expected);
    end
end

function A = load_auxiliary_transect_columns( ...
        inputFolder, caseBaseName, csvSuffix, profileIndex, rowCoord, ...
        transectN_m, tolerance)
    sourceFileName = string(caseBaseName) + string(csvSuffix);
    sourceFile = fullfile(inputFolder, sourceFileName);
    if ~isfile(sourceFile)
        error(['Required auxiliary transect CSV was not found for %s: %s. ', ...
            'Run calling_analysis_velocity_export.m for this case first.'], ...
            string(caseBaseName), sourceFile);
    end

    source = readtable(sourceFile, ...
        'TextType', 'string', ...
        'VariableNamingRule', 'preserve');
    requiredColumns = [ ...
        "profileIndex", "rowCoord", "rowCell", "transectN_m", ...
        "nearestRowIndex", "crossSectionX_m", "crossSectionY_m", ...
        "crossSectionBedElevation_m", ...
        "crossSectionInitialMedianWSE_m", "crossSectionDepth_m", ...
        "crossSectionDepthStatus", "nIntersectedCells", ...
        "nFiniteVelocityCells"];
    missingColumns = setdiff(requiredColumns, ...
        string(source.Properties.VariableNames));
    if ~isempty(missingColumns)
        error('Auxiliary transect CSV %s is missing column(s): %s', ...
            sourceFile, char(strjoin(missingColumns, ', ')));
    end

    sourceProfileIndex = round(double(source.profileIndex));
    if any(~isfinite(sourceProfileIndex)) || ...
            numel(unique(sourceProfileIndex)) ~= numel(sourceProfileIndex)
        error('Auxiliary transect CSV has invalid or duplicate profileIndex values: %s', ...
            sourceFile);
    end

    [isPresent, sourceRow] = ismember(profileIndex(:), sourceProfileIndex);
    if ~all(isPresent)
        missingProfileIndex = profileIndex(~isPresent);
        error('Auxiliary transect CSV %s is missing profileIndex value(s): %s', ...
            sourceFile, char(strjoin(string(missingProfileIndex), ', ')));
    end
    source = source(sourceRow, :);

    validate_numeric_vector(double(source.rowCoord), rowCoord, tolerance, ...
        sourceFileName + " rowCoord");
    validate_numeric_vector(double(source.transectN_m), transectN_m, ...
        tolerance, sourceFileName + " transectN_m");

    A = struct;
    A.rowCell = double(source.rowCell);
    A.nearestRowIndex = double(source.nearestRowIndex);
    A.crossSectionX_m = double(source.crossSectionX_m);
    A.crossSectionY_m = double(source.crossSectionY_m);
    A.crossSectionBedElevation_m = ...
        double(source.crossSectionBedElevation_m);
    A.crossSectionInitialMedianWSE_m = ...
        double(source.crossSectionInitialMedianWSE_m);
    A.crossSectionDepth_m = double(source.crossSectionDepth_m);
    A.crossSectionDepthStatus = string(source.crossSectionDepthStatus);
    A.nIntersectedCells = double(source.nIntersectedCells);
    A.nFiniteVelocityCells = double(source.nFiniteVelocityCells);
end

function validate_numeric_vector(actual, expected, tolerance, label)
    actual = double(actual(:));
    expected = double(expected(:));
    if numel(actual) ~= numel(expected)
        error('%s length mismatch: auxiliary %d, selected profiles %d.', ...
            label, numel(actual), numel(expected));
    end

    sameNaN = isnan(actual) & isnan(expected);
    bothFinite = isfinite(actual) & isfinite(expected);
    scale = max([ones(numel(actual), 1), abs(actual), abs(expected)], ...
        [], 2);
    closeFinite = bothFinite & ...
        abs(actual - expected) <= tolerance .* scale;
    if ~all(sameNaN | closeFinite)
        firstMismatch = find(~(sameNaN | closeFinite), 1, 'first');
        error('%s mismatch at joined row %d: auxiliary %.17g, selected %.17g.', ...
            label, firstMismatch, actual(firstMismatch), expected(firstMismatch));
    end
end

function text = ternary_text(condition, trueText, falseText)
    if condition
        text = trueText;
    else
        text = falseText;
    end
end
