function outputs = generate_case_transect_analysis(archiveRoot, baseName, outputFolder)
%GENERATE_CASE_TRANSECT_ANALYSIS Rebuild one real-case transect product.
%   OUTPUTS = GENERATE_CASE_TRANSECT_ANALYSIS(ARCHIVEROOT, BASENAME,
%   OUTPUTFOLDER) reads the deposited checkpoint, velocity input and
%   cross-section files and recreates the auxiliary per-transect CSV used
%   by the selected-profile exporter.
arguments
    archiveRoot (1,1) string
    baseName (1,1) string
    outputFolder (1,1) string
end

% Focused first-10 m transect analysis using observed 0.5 s tracked velocity for inverse wavelength-depth estimates.
%
% This script deliberately removes the older central-transect-only,
% instantaneous-velocity, and plotting/export sections that are not required
% for this analysis.
%
% Workflow:
%   1) Run the refined checkpoint selector over the first 10 m from the first
%      WSE-populated row.
%   2) Generate velocity maps using ONLY the 0.5 s vector-tracking input file.
%   3) For every selector transect across the first 10 m:
%        a) extract the median 0.5 s tracked streamwise velocity from cells
%           intersected by that transect;
%        b) get wavelength from the accepted/best WSE map using autocorrelation and PMUSIC;
%        c) get wavelength from the initial residual map using autocorrelation and PMUSIC;
%        d) estimate depth from observed dt=0.5 s median tracked velocity
%           plus each wavelength source, including PMUSIC wavelengths;
%        e) estimate deep-water gravity-wave velocity from each autocorrelation
%           and PMUSIC wavelength;
%        f) retain WSE-wavelength, initial-residual-wavelength, WSE-PMUSIC,
%           and initial-residual-PMUSIC depth estimates.
%
% PMUSIC wavelengths are exported as additional wavelength estimates and are
% also used to calculate additional inverse-depth estimates using the same
% observed-velocity wavelength-depth relation as the autocorrelation wavelengths.
%
% Main outputs:
%   <baseName>_first10m_transect_depth_from_observed_velocity_with_deep_water_velocity.mat
%   <baseName>_first10m_transect_depth_from_observed_velocity_with_deep_water_velocity.csv
%
% The output table also includes cross-section-derived bed elevation and
% water depth at each transect location, using crossSectionCsvFile.
% Per-transect wave amplitudes are also exported:
%   acceptedWseAmplitude_m = 0.5*(P98-P2) of detrended accepted WSE profile
%   initialResidualAmplitude = 0.5*(P98-P2) of detrended initial residual profile
%   wsePmusicWavelength_m and initialResidualPmusicWavelength_m are PMUSIC-based
%   wavelength estimates from the same transect signals.
%   depthFromWsePmusicWavelength_m and depthFromInitialResidualPmusicWavelength_m
%   use observedVelocityTrackedMedian_mps with the PMUSIC wavelength estimates.

thisScriptDir = fileparts(mfilename('fullpath'));
if ~isempty(thisScriptDir)
    addpath(thisScriptDir, '-begin');
end
clear KLT_select_wave_checkpoint_refined KLT_select_wave_checkpoint_refined_first10m
rehash;

fprintf('\nUsing selector: %s\n\n', which('KLT_select_wave_checkpoint_refined'));

%% ------------------------------------------------------------------------
% User settings
%% ------------------------------------------------------------------------
archiveRoot = string(java.io.File(archiveRoot).getCanonicalPath());
inputFolder = fullfile(archiveRoot, "videos", "inputs");
checkpointFolder = fullfile(archiveRoot, "videos", "checkpoints");
if ~isfolder(inputFolder) || ~isfolder(checkpointFolder)
    error('Expected videos/inputs and videos/checkpoints below %s.', archiveRoot);
end
if ~isfolder(outputFolder)
    mkdir(outputFolder);
end

checkpointFile = fullfile(checkpointFolder, baseName + "_checkpoint.mat");
velocityInputsFile = fullfile(inputFolder, baseName + "_velocity_inputs.mat");
if ~isfile(checkpointFile) || ~isfile(velocityInputsFile)
    error('Missing checkpoint or velocity-input file for %s.', baseName);
end
velocityInputsDt_s  = 0.5;

outputMatFile = fullfile(outputFolder, baseName + ...
    "_first10m_transect_depth_from_observed_velocity_with_deep_water_velocity_pmusic_depths.mat");
outputCsvFile = fullfile(outputFolder, baseName + ...
    "_first10m_transect_depth_from_observed_velocity_with_deep_water_velocity_pmusic_depths.csv");

saveOutputs = true;
saveMatOutput = false;


% Cross-section file used for independent bed-elevation / water-depth extraction.
% The x/y coordinates are projected into the same solver s/n coordinate system
% as the WSE grid, then interpolated at each first-10 m transect location.
crossSectionCsvFile = fullfile(inputFolder, 'cross_section.csv');
crossSectionCsvFallbackFiles = {char(crossSectionCsvFile)};
adjustCrossSectionToWseExtent = true;
crossSectionAlignmentDirection = 'solver_n';
waterDepthWseMapAA = 1;
waterDepthWseMapIndex = 1;

% First-10 m selector settings.
% These are matched to calling_syn_analysis so the same checkpoint should
% produce the same selected stop map.  In particular, do NOT force
% minPopulatedCellsPerRow = 1 here unless it is also used in the comparison
% script, because that can move the first-10 m analysis band.
selectorOpts = struct;
selectorOpts.stabilityAssessmentMode = 'first_populated_10m_from_camera';
selectorOpts.firstSensedLength_m = 10;
selectorOpts.cameraRowOrder = 'ascending';
selectorOpts.displaySummary = true;

% PMUSIC configuration, matched to calling_syn_analysis.  PMUSIC is reported
% by the selector, but the stop-map decision is primarily controlled by
% amplitude/wavelength stability, shape change, residual metrics and the
% amplitude-plateau checks below.
selectorOpts.usePmusicWavelength = true;
selectorOpts.pmusicModelOrder = 2;
selectorOpts.pmusicNfft = 2048;
selectorOpts.pmusicWindowLength = [];
selectorOpts.pmusicOverlap = [];

% Robust stop-map selection, matched to calling_syn_analysis.
selectorOpts.requireForwardAmplitudePlateau = true;
selectorOpts.weakForwardWindow = 5;
selectorOpts.weakMaxForwardAmpGrowth = 0.05;
selectorOpts.weakMaxForwardAmpAbsGrowth_m = 0.003;
selectorOpts.strongForwardWindow = 5;
selectorOpts.strongMaxForwardAmpGrowth = 0.05;
selectorOpts.strongMaxForwardAmpAbsGrowth_m = 0.003;
selectorOpts.noPlateauFallbackMode = 'latest_supported';

% Leave this unset so the selector default is used, matching calling_syn_analysis.
% selectorOpts.minPopulatedCellsPerRow = 1;

% 0.5 s velocity-map settings.
trackedVelOpts = struct;
trackedVelOpts.pathPopulation = 'all';
trackedVelOpts.binAnchor = 'A';
trackedVelOpts.dt_s = velocityInputsDt_s;
trackedVelOpts.keepProjectedPaths = false;
trackedVelOpts.saveMatFile = '';

% Wavelength/amplitude/depth analysis settings.
% Inverse depth estimates use observedVelocityTrackedMedian_mps plus each
% wavelength source: WSE wavelength, initial-residual wavelength, WSE-PMUSIC wavelength,
% and initial-residual-PMUSIC wavelength.
analysisOpts = struct;
analysisOpts.depthGrid_m = 0.01:0.01:10;
analysisOpts.reasonableDepthRange_m = [0.05 5.00];
analysisOpts.maxWavelengthMismatchFraction = 0.10;
analysisOpts.minWavelength_m = 0.30;
analysisOpts.maxWavelength_m = Inf;
analysisOpts.maxLagFraction = 0.75;
analysisOpts.minPeakCorrelation = 0.10;
analysisOpts.detrendOrder = 1;
analysisOpts.ampLowPct = 2;
analysisOpts.ampHighPct = 98;
analysisOpts.usePmusicWavelength = true;
analysisOpts.pmusicModelOrder = 2;          % 2 for one real sinusoidal component
analysisOpts.pmusicNfft = 2048;
analysisOpts.pmusicWindowLength = [];       % [] chooses a data-dependent value
analysisOpts.pmusicOverlap = [];            % [] chooses 50% overlap

%% ------------------------------------------------------------------------
% Load checkpoint and run first-10 m WSE selector
%% ------------------------------------------------------------------------
S = load(checkpointFile, 'checkpoint');
if ~isfield(S, 'checkpoint')
    error('checkpointFile must contain a variable named checkpoint.');
end
checkpoint = S.checkpoint;

if isfield(checkpoint, 'aa') && isfinite(checkpoint.aa)
    aa = checkpoint.aa;
else
    aa = 1;
end

% Use the refined selector explicitly rather than falling back to a different
% first10m selector version.  This removes path/order ambiguity and makes the
% stop-map selection directly comparable with calling_syn_analysis.
selectorFunction = 'KLT_select_wave_checkpoint_refined';
if exist(selectorFunction, 'file') ~= 2
    error('KLT_select_wave_checkpoint_refined.m was not found on the MATLAB path.');
end

[result, mapTable, profileTable] = KLT_select_wave_checkpoint_refined(checkpoint, selectorOpts);

fprintf('\n=== Selector diagnostic ===\n');
fprintf('Selector function: %s\n', which('KLT_select_wave_checkpoint_refined'));
fprintf('Selected map index: %d\n', result.selectedMapIndex);
fprintf('Stop method: %s\n', char(string(result.stopMethod)));
if isfield(result, 'firstPopulatedRowIndex')
    fprintf('First populated row index: %d\n', result.firstPopulatedRowIndex);
end
if isfield(result, 'stabilityBandRows') && ~isempty(result.stabilityBandRows)
    fprintf('Stability band rows: %d to %d, n = %d\n', ...
        min(result.stabilityBandRows), max(result.stabilityBandRows), numel(result.stabilityBandRows));
end
if isfield(result, 'nAmplitudeProfiles')
    fprintf('nAmplitudeProfiles: %d\n', result.nAmplitudeProfiles);
end
if isfield(result, 'nWavelengthProfiles')
    fprintf('nWavelengthProfiles: %d\n', result.nWavelengthProfiles);
end
if isfield(result, 'stopInfo')
    disp(result.stopInfo);
end

acceptedMapIndex = resolve_accepted_map_index(checkpoint, aa, result, mapTable);
acceptedWseMap = double(checkpoint.wse_map{aa, acceptedMapIndex});
initialWseMap = double(checkpoint.wse_map{aa, 1});

[initialResidualMap, residualInfo] = get_initial_residual_map(checkpoint, aa, size(acceptedWseMap));

%% ------------------------------------------------------------------------
% Generate velocity maps using ONLY the 0.5 s vector-tracking inputs
%% ------------------------------------------------------------------------
velocityFunction = choose_existing_function({ ...
    'make_velocity_maps_from_checkpoint', ...
    'KLT_make_velocity_maps_from_checkpoint'});

trackedVelOpts.acceptedMapIndex = acceptedMapIndex;

velocityOutTracked = feval(velocityFunction, ...
    checkpoint, ...
    velocityInputsFile, ...
    result, ...
    mapTable, ...
    trackedVelOpts);

initialVelocityMapTracked_mps  = velocityOutTracked.start.u_streamwise_mps;
acceptedVelocityMapTracked_mps = velocityOutTracked.accepted.u_streamwise_mps;
velocityDifferenceMapTracked_mps = acceptedVelocityMapTracked_mps - initialVelocityMapTracked_mps;

%% ------------------------------------------------------------------------
% Analyse every first-10 m transect
%% ------------------------------------------------------------------------
if ~isfield(result, 'rowCoord') || ~isfield(result, 'colCoord') || ~isfield(result, 'distance_m')
    error('Selector result must contain rowCoord, colCoord, and distance_m.');
end

rowCoord = result.rowCoord(:);
colCoord = result.colCoord(:);
distance_m = result.distance_m(:);
nProfiles = numel(rowCoord);

% Convert each selector transect from fractional row index into the solver
% cross-stream physical coordinate n.  These n-locations are the exact
% locations used for interpolating the cross-section bed/depth values.
yiForTransects = checkpoint.yi(:);
transectN_m = interp1((1:numel(yiForTransects)).', yiForTransects, rowCoord, 'linear', NaN);
nearestRowIndex = nan(nProfiles, 1);
for pp = 1:nProfiles
    [~, nearestRowIndex(pp)] = min(abs(yiForTransects - transectN_m(pp)));
end

[transectCrossSectionDepthTable, crossSectionOverlay, crossSectionDepthInfo] = ...
    compute_cross_section_depths_for_transects( ...
        checkpoint, transectN_m, rowCoord, nearestRowIndex, ...
        crossSectionCsvFile, crossSectionCsvFallbackFiles, ...
        adjustCrossSectionToWseExtent, crossSectionAlignmentDirection, ...
        waterDepthWseMapAA, waterDepthWseMapIndex);

transectAnalysis = repmat(empty_transect_struct(), nProfiles, 1);

for pp = 1:nProfiles
    cells = get_transect_cells(rowCoord, colCoord, distance_m, pp, size(acceptedWseMap));

    initialVelocityVals = initialVelocityMapTracked_mps(cells.linearIndices);
    velocityVals = acceptedVelocityMapTracked_mps(cells.linearIndices);
    wseVals = acceptedWseMap(cells.linearIndices);
    finiteInitialVel = isfinite(initialVelocityVals) & isfinite(wseVals);
    finiteVel = isfinite(velocityVals) & isfinite(wseVals);

    if any(finiteInitialVel)
        initialVelocity_mps = median( ...
            initialVelocityVals(finiteInitialVel), 'omitnan');
    else
        initialVelocity_mps = NaN;
    end

    if any(finiteVel)
        observedVelocity_mps = median(velocityVals(finiteVel), 'omitnan');
    else
        observedVelocity_mps = NaN;
    end

    [wseWavelength_m, wseWavelengthSource, wsePeakR] = get_wse_wavelength_for_profile( ...
        profileTable, acceptedMapIndex, pp, acceptedWseMap, rowCoord, colCoord, distance_m, analysisOpts);

    [acceptedWseAmplitude_m, acceptedWseAmplitudeSource, acceptedWseAmplitudeStatus] = get_wse_amplitude_for_profile( ...
        profileTable, acceptedMapIndex, pp, acceptedWseMap, rowCoord, colCoord, distance_m, analysisOpts);

    residualW = estimate_map_wavelength_along_profile( ...
        initialResidualMap, rowCoord, colCoord, distance_m, pp, analysisOpts);
    residualWavelength_m = residualW.wavelength_m;

    wsePmusicW = estimate_map_wavelength_pmusic_along_profile( ...
        acceptedWseMap, rowCoord, colCoord, distance_m, pp, analysisOpts);
    residualPmusicW = estimate_map_wavelength_pmusic_along_profile( ...
        initialResidualMap, rowCoord, colCoord, distance_m, pp, analysisOpts);

    residualA = estimate_map_amplitude_along_profile( ...
        initialResidualMap, rowCoord, colCoord, distance_m, pp, analysisOpts);

    vDeepWaterFromWse = velocity_from_deep_water_gravity_wave(wseWavelength_m);
    vDeepWaterFromResidual = velocity_from_deep_water_gravity_wave(residualWavelength_m);
    vDeepWaterFromWsePmusic = velocity_from_deep_water_gravity_wave(wsePmusicW.wavelength_m);
    vDeepWaterFromResidualPmusic = velocity_from_deep_water_gravity_wave(residualPmusicW.wavelength_m);

    depthFromWse = estimate_depth_from_known_wavelength_filtered( ...
        observedVelocity_mps, wseWavelength_m, analysisOpts);
    depthFromResidual = estimate_depth_from_known_wavelength_filtered( ...
        observedVelocity_mps, residualWavelength_m, analysisOpts);
    depthFromWsePmusic = estimate_depth_from_known_wavelength_filtered( ...
        observedVelocity_mps, wsePmusicW.wavelength_m, analysisOpts);
    depthFromResidualPmusic = estimate_depth_from_known_wavelength_filtered( ...
        observedVelocity_mps, residualPmusicW.wavelength_m, analysisOpts);

    transectAnalysis(pp).profileIndex = pp;
    transectAnalysis(pp).rowCoord = rowCoord(pp);
    transectAnalysis(pp).rowCell = cells.rowCell;
    transectAnalysis(pp).transectN_m = transectN_m(pp);
    transectAnalysis(pp).nearestRowIndex = nearestRowIndex(pp);
    transectAnalysis(pp).crossSectionX_m = transectCrossSectionDepthTable.crossSectionX_m(pp);
    transectAnalysis(pp).crossSectionY_m = transectCrossSectionDepthTable.crossSectionY_m(pp);
    transectAnalysis(pp).crossSectionBedElevation_m = transectCrossSectionDepthTable.crossSectionBedElevation_m(pp);
    transectAnalysis(pp).crossSectionInitialMedianWSE_m = transectCrossSectionDepthTable.initialMedianWSE_m(pp);
    transectAnalysis(pp).crossSectionDepth_m = transectCrossSectionDepthTable.crossSectionDepth_m(pp);
    transectAnalysis(pp).crossSectionDepthStatus = char(transectCrossSectionDepthTable.crossSectionDepthStatus(pp));
    transectAnalysis(pp).nIntersectedCells = numel(cells.linearIndices);
    transectAnalysis(pp).nFiniteVelocityCells = nnz(finiteVel);
    transectAnalysis(pp).initialVelocityTrackedMedian_mps = ...
        initialVelocity_mps;
    transectAnalysis(pp).acceptedVelocityTrackedMedian_mps = ...
        observedVelocity_mps;
    transectAnalysis(pp).observedVelocityTrackedMedian_mps = observedVelocity_mps;

    transectAnalysis(pp).wseWavelength_m = wseWavelength_m;
    transectAnalysis(pp).acceptedWseAmplitude_m = acceptedWseAmplitude_m;
    transectAnalysis(pp).acceptedWseAmplitudeStatus = acceptedWseAmplitudeStatus;
    transectAnalysis(pp).initialResidualWavelength_m = residualWavelength_m;
    transectAnalysis(pp).initialResidualWavelengthStatus = residualW.status;
    transectAnalysis(pp).wsePmusicWavelength_m = wsePmusicW.wavelength_m;
    transectAnalysis(pp).wsePmusicWavelengthStatus = wsePmusicW.status;
    transectAnalysis(pp).wsePmusicPeakFrequency_cpm = wsePmusicW.peakFrequency_cpm;
    transectAnalysis(pp).wsePmusicPeakPower = wsePmusicW.peakPower;
    transectAnalysis(pp).initialResidualPmusicWavelength_m = residualPmusicW.wavelength_m;
    transectAnalysis(pp).initialResidualPmusicWavelengthStatus = residualPmusicW.status;
    transectAnalysis(pp).initialResidualPmusicPeakFrequency_cpm = residualPmusicW.peakFrequency_cpm;
    transectAnalysis(pp).initialResidualPmusicPeakPower = residualPmusicW.peakPower;
    transectAnalysis(pp).initialResidualAmplitude = residualA.amplitude;
    transectAnalysis(pp).initialResidualAmplitudeUnits = residualA.units;
    transectAnalysis(pp).initialResidualAmplitudeStatus = residualA.status;

    transectAnalysis(pp).velocityFromWseWavelengthDeepWaterGravity_mps = vDeepWaterFromWse.velocity_mps;
    transectAnalysis(pp).velocityFromWseWavelengthDeepWaterGravityStatus = vDeepWaterFromWse.status;
    transectAnalysis(pp).velocityFromResidualWavelengthDeepWaterGravity_mps = vDeepWaterFromResidual.velocity_mps;
    transectAnalysis(pp).velocityFromResidualWavelengthDeepWaterGravityStatus = vDeepWaterFromResidual.status;
    transectAnalysis(pp).velocityFromWsePmusicWavelengthDeepWaterGravity_mps = vDeepWaterFromWsePmusic.velocity_mps;
    transectAnalysis(pp).velocityFromWsePmusicWavelengthDeepWaterGravityStatus = vDeepWaterFromWsePmusic.status;
    transectAnalysis(pp).velocityFromResidualPmusicWavelengthDeepWaterGravity_mps = vDeepWaterFromResidualPmusic.velocity_mps;
    transectAnalysis(pp).velocityFromResidualPmusicWavelengthDeepWaterGravityStatus = vDeepWaterFromResidualPmusic.status;

    transectAnalysis(pp).depthFromWseWavelength_m = depthFromWse.depthEstimate_m;
    transectAnalysis(pp).depthFromWseWavelengthStatus = depthFromWse.status;
    transectAnalysis(pp).depthFromInitialResidualWavelength_m = depthFromResidual.depthEstimate_m;
    transectAnalysis(pp).depthFromInitialResidualWavelengthStatus = depthFromResidual.status;
    transectAnalysis(pp).depthFromWsePmusicWavelength_m = depthFromWsePmusic.depthEstimate_m;
    transectAnalysis(pp).depthFromWsePmusicWavelengthStatus = depthFromWsePmusic.status;
    transectAnalysis(pp).depthFromInitialResidualPmusicWavelength_m = depthFromResidualPmusic.depthEstimate_m;
    transectAnalysis(pp).depthFromInitialResidualPmusicWavelengthStatus = depthFromResidualPmusic.status;

    transectAnalysis(pp).cellRows = cells.rows;
    transectAnalysis(pp).cellCols = cells.cols;
    transectAnalysis(pp).cellLinearIndices = cells.linearIndices;
    transectAnalysis(pp).velocityCellsTracked_mps = velocityVals(:);
end

transectAnalysisTable = struct2table(remove_large_fields_for_table(transectAnalysis));
% The table is intentionally concise. Depth estimates are calculated from
% observedVelocityTrackedMedian_mps plus each wavelength source.
% Raw transect cell arrays remain in transectAnalysis inside the MAT file only.

%% ------------------------------------------------------------------------
% Coordinates and summary
%% ------------------------------------------------------------------------
xCoordinates_m = checkpoint.xi(:).';
yCoordinates_m = checkpoint.yi(:);

if isfield(checkpoint, 'X_rot') && isfield(checkpoint, 'Y_rot') && ...
        isequal(size(checkpoint.X_rot), size(acceptedWseMap)) && ...
        isequal(size(checkpoint.Y_rot), size(acceptedWseMap))
    XCoordinates_m = checkpoint.X_rot;
    YCoordinates_m = checkpoint.Y_rot;
else
    [XCoordinates_m, YCoordinates_m] = meshgrid(xCoordinates_m, yCoordinates_m);
end

summary = struct;
summary.nTransects = nProfiles;
summary.nFiniteObservedVelocity = nnz(isfinite(transectAnalysisTable.observedVelocityTrackedMedian_mps));
summary.nFiniteCrossSectionDepths = nnz(isfinite(transectAnalysisTable.crossSectionDepth_m));
summary.nValidEstimatedDepthsFromWseWavelength = nnz(isfinite(transectAnalysisTable.depthFromWseWavelength_m));
summary.nValidEstimatedDepthsFromInitialResidualWavelength = nnz(isfinite(transectAnalysisTable.depthFromInitialResidualWavelength_m));
summary.nValidEstimatedDepthsFromWsePmusicWavelength = nnz(isfinite(transectAnalysisTable.depthFromWsePmusicWavelength_m));
summary.nValidEstimatedDepthsFromInitialResidualPmusicWavelength = nnz(isfinite(transectAnalysisTable.depthFromInitialResidualPmusicWavelength_m));
summary.nValidDeepWaterGravityWseVelocity = nnz(isfinite(transectAnalysisTable.velocityFromWseWavelengthDeepWaterGravity_mps));
summary.nValidDeepWaterGravityInitialResidualVelocity = nnz(isfinite(transectAnalysisTable.velocityFromResidualWavelengthDeepWaterGravity_mps));
summary.nValidDeepWaterGravityWsePmusicVelocity = nnz(isfinite(transectAnalysisTable.velocityFromWsePmusicWavelengthDeepWaterGravity_mps));
summary.nValidDeepWaterGravityInitialResidualPmusicVelocity = nnz(isfinite(transectAnalysisTable.velocityFromResidualPmusicWavelengthDeepWaterGravity_mps));
summary.nValidPmusicWseWavelengths = nnz(isfinite(transectAnalysisTable.wsePmusicWavelength_m));
summary.nValidPmusicInitialResidualWavelengths = nnz(isfinite(transectAnalysisTable.initialResidualPmusicWavelength_m));
summary.medianObservedTrackedVelocity_mps = median(transectAnalysisTable.observedVelocityTrackedMedian_mps, 'omitnan');
summary.medianWseWavelength_m = median(transectAnalysisTable.wseWavelength_m, 'omitnan');
summary.medianAcceptedWseAmplitude_m = median(transectAnalysisTable.acceptedWseAmplitude_m, 'omitnan');
summary.medianInitialResidualWavelength_m = median(transectAnalysisTable.initialResidualWavelength_m, 'omitnan');
summary.medianWsePmusicWavelength_m = median(transectAnalysisTable.wsePmusicWavelength_m, 'omitnan');
summary.medianInitialResidualPmusicWavelength_m = median(transectAnalysisTable.initialResidualPmusicWavelength_m, 'omitnan');
summary.medianInitialResidualAmplitude = median(transectAnalysisTable.initialResidualAmplitude, 'omitnan');
summary.medianCrossSectionDepth_m = median(transectAnalysisTable.crossSectionDepth_m, 'omitnan');
summary.medianDepthFromWseWavelength_m = median(transectAnalysisTable.depthFromWseWavelength_m, 'omitnan');
summary.medianDepthFromInitialResidualWavelength_m = median(transectAnalysisTable.depthFromInitialResidualWavelength_m, 'omitnan');
summary.medianDepthFromWsePmusicWavelength_m = median(transectAnalysisTable.depthFromWsePmusicWavelength_m, 'omitnan');
summary.medianDepthFromInitialResidualPmusicWavelength_m = median(transectAnalysisTable.depthFromInitialResidualPmusicWavelength_m, 'omitnan');
summary.medianVelocityFromWseWavelengthDeepWaterGravity_mps = median(transectAnalysisTable.velocityFromWseWavelengthDeepWaterGravity_mps, 'omitnan');
summary.medianVelocityFromResidualWavelengthDeepWaterGravity_mps = median(transectAnalysisTable.velocityFromResidualWavelengthDeepWaterGravity_mps, 'omitnan');
summary.medianVelocityFromWsePmusicWavelengthDeepWaterGravity_mps = median(transectAnalysisTable.velocityFromWsePmusicWavelengthDeepWaterGravity_mps, 'omitnan');
summary.medianVelocityFromResidualPmusicWavelengthDeepWaterGravity_mps = median(transectAnalysisTable.velocityFromResidualPmusicWavelengthDeepWaterGravity_mps, 'omitnan');

exportInfo = struct;
exportInfo.baseName = baseName;
exportInfo.checkpointFile = checkpointFile;
exportInfo.velocityInputsFile = velocityInputsFile;
exportInfo.velocityInputsDt_s = velocityInputsDt_s;
exportInfo.crossSectionCsvFile = crossSectionCsvFile;
exportInfo.crossSectionFileUsed = crossSectionDepthInfo.fileUsed;
exportInfo.crossSectionAdjustedToWseExtent = crossSectionDepthInfo.adjustedToWseExtent;
exportInfo.crossSectionAlignmentDirection = crossSectionDepthInfo.alignmentDirection;
exportInfo.waterDepthWseMapAA = waterDepthWseMapAA;
exportInfo.waterDepthWseMapIndex = waterDepthWseMapIndex;
exportInfo.selectorFunction = selectorFunction;
exportInfo.velocityFunction = velocityFunction;
exportInfo.aa = aa;
exportInfo.initialWseMapIndex = 1;
exportInfo.acceptedMapIndex = acceptedMapIndex;
exportInfo.initialResidualMapSource = residualInfo.source;
exportInfo.initialResidualMapIndex = residualInfo.index;
exportInfo.velocitySource = '0.5 s vector-tracking inputs only';
exportInfo.transectSet = 'all selector transects across first 10 m';
exportInfo.tableOutputMode = 'inverse depth estimates from observed dt=0.5 s velocity and WSE/residual wavelengths';
exportInfo.observedVelocityDefinition = 'median accepted-WSE-projected 0.5 s tracked streamwise velocity in cells intersected by each transect';
exportInfo.wseWavelengthDefinition = 'accepted/best WSE-map autocorrelation wavelength per transect';
exportInfo.acceptedWseAmplitudeDefinition = 'accepted/best WSE-map transect amplitude: 0.5*(P98-P2) after linear detrending; sourced from profileTable.amplitude98_2_m when available';
exportInfo.initialResidualWavelengthDefinition = 'initial residual map sampled along the same transect; first positive spatial-autocorrelation peak';
exportInfo.pmusicWavelengthDefinition = ['PMUSIC wavelength estimate from the same detrended transect signal. ', ...
    'The PMUSIC peak frequency is converted to wavelength as lambda = 1/f, where f is cycles per metre. ', ...
    'These PMUSIC wavelengths are also used for additional inverse-depth and deep-water velocity estimates.'];
exportInfo.initialResidualAmplitudeDefinition = 'initial residual map transect amplitude: 0.5*(P98-P2) after linear detrending; units follow the residual map source';
exportInfo.crossSectionDepthDefinition = 'known water depth from crossSectionCsvFile interpolated at each transect n-location: initial median WSE minus cross-section bed elevation';
exportInfo.estimatedDepthDefinition = ['inverse depth from observed dt=0.5 s median velocity and measured wavelength ', ...
    'using the supplied wavelength-depth relation; calculated separately for autocorrelation and PMUSIC wavelengths'];
exportInfo.deepWaterGravityVelocityDefinition = 'velocity calculated from autocorrelation and PMUSIC wavelengths using U_s = sqrt(g*lambda/(2*pi))';
exportInfo.estimatedDepthDiscardRule = sprintf(['discard estimated depth if wavelength is outside the finite curve range, ', ...
    'if mismatch > %.3g, or if raw depth is outside [%.3g %.3g] m'], ...
    analysisOpts.maxWavelengthMismatchFraction, analysisOpts.reasonableDepthRange_m(1), analysisOpts.reasonableDepthRange_m(2));
exportInfo.saveOutputs = saveOutputs;
exportInfo.velocityUnits = 'm/s';
exportInfo.coordinateDescription = ['xCoordinates_m and yCoordinates_m are the 1-D checkpoint grid axes; ', ...
    'XCoordinates_m and YCoordinates_m are Ny-by-Nx map-coordinate arrays matching the exported maps.'];

%% ------------------------------------------------------------------------
% Save outputs
%% ------------------------------------------------------------------------
if saveOutputs
    writetable(transectAnalysisTable, outputCsvFile);

    if saveMatOutput
        save(outputMatFile, ...
        'baseName', ...
        'initialResidualMap', ...
        'initialWseMap', ...
        'acceptedWseMap', ...
        'initialVelocityMapTracked_mps', ...
        'acceptedVelocityMapTracked_mps', ...
        'velocityDifferenceMapTracked_mps', ...
        'transectAnalysis', ...
        'transectAnalysisTable', ...
        'transectCrossSectionDepthTable', ...
        'crossSectionOverlay', ...
        'crossSectionDepthInfo', ...
        'summary', ...
        'xCoordinates_m', ...
        'yCoordinates_m', ...
        'XCoordinates_m', ...
        'YCoordinates_m', ...
        'result', ...
        'mapTable', ...
        'profileTable', ...
        'velocityOutTracked', ...
        'analysisOpts', ...
        'exportInfo', ...
            '-v7.3');
    end

    fprintf('Saved first-10 m transect analysis CSV file to %s\n', outputCsvFile);
else
    fprintf('\nSaving disabled: saveOutputs = false. No MAT/CSV files written.\n');
end
fprintf('Accepted WSE map index: %d\n', acceptedMapIndex);
fprintf('Initial residual source: %s, index: %d\n', residualInfo.source, residualInfo.index);
fprintf('0.5 s velocity input file: %s\n', velocityInputsFile);
fprintf('Transects analysed: %d\n', summary.nTransects);
fprintf('Finite observed 0.5 s velocities: %d\n', summary.nFiniteObservedVelocity);
fprintf('Finite cross-section depths at transect locations: %d\n', summary.nFiniteCrossSectionDepths);
fprintf('Valid PMUSIC WSE wavelengths: %d\n', summary.nValidPmusicWseWavelengths);
fprintf('Valid PMUSIC initial-residual wavelengths: %d\n', summary.nValidPmusicInitialResidualWavelengths);
fprintf('Valid PMUSIC WSE wavelength depth estimates: %d\n', summary.nValidEstimatedDepthsFromWsePmusicWavelength);
fprintf('Valid PMUSIC initial-residual wavelength depth estimates: %d\n', summary.nValidEstimatedDepthsFromInitialResidualPmusicWavelength);
fprintf('Median accepted WSE amplitude: %.5g m\n', summary.medianAcceptedWseAmplitude_m);
fprintf('Median initial residual amplitude: %.5g\n', summary.medianInitialResidualAmplitude);
fprintf('Cross-section file used: %s\n', crossSectionDepthInfo.fileUsed);

outputs = struct;
outputs.csvFile = char(outputCsvFile);
outputs.matFile = '';
outputs.transectAnalysisTable = transectAnalysisTable;
outputs.result = result;
outputs.mapTable = mapTable;
outputs.profileTable = profileTable;
outputs.summary = summary;
end

%% ========================================================================
% Local helper functions
%% ========================================================================

function [T, C, info] = compute_cross_section_depths_for_transects( ...
    checkpoint, transectN_m, rowCoord, nearestRowIndex, ...
    crossSectionCsvFile, crossSectionCsvFallbackFiles, ...
    adjustToWseExtent, alignmentDirection, waterDepthWseMapAA, waterDepthWseMapIndex)
% Compute cross-section-derived depth at exactly the same n-locations as
% the analysed first-10 m transects.

    nProfiles = numel(transectN_m);
    C = build_cross_section_overlay_for_depth( ...
        crossSectionCsvFile, crossSectionCsvFallbackFiles, ...
        checkpoint, adjustToWseExtent, alignmentDirection);

    info = struct;
    info.fileRequested = crossSectionCsvFile;
    info.fileUsed = '';
    info.loaded = false;
    info.adjustedToWseExtent = false;
    info.alignmentDirection = char(alignmentDirection);
    info.status = 'not evaluated';

    profileIndex = (1:nProfiles).';
    crossSectionX_m = nan(nProfiles, 1);
    crossSectionY_m = nan(nProfiles, 1);
    crossSectionBedElevation_m = nan(nProfiles, 1);
    initialMedianWSE_m = nan(nProfiles, 1);
    crossSectionDepth_m = nan(nProfiles, 1);
    status = repmat("not evaluated", nProfiles, 1);

    if isfield(C, 'fileUsed'), info.fileUsed = C.fileUsed; end
    if isfield(C, 'loaded'), info.loaded = C.loaded; end
    if isfield(C, 'adjustedToWseExtent'), info.adjustedToWseExtent = C.adjustedToWseExtent; end
    if isfield(C, 'adjustmentDirectionName'), info.alignmentDirection = C.adjustmentDirectionName; end

    if ~isfield(C, 'loaded') || ~C.loaded
        status(:) = "cross-section not loaded";
        info.status = 'cross-section not loaded';
        T = make_cross_section_depth_table(profileIndex, transectN_m, rowCoord, nearestRowIndex, ...
            crossSectionX_m, crossSectionY_m, crossSectionBedElevation_m, initialMedianWSE_m, crossSectionDepth_m, status);
        return
    end

    if ~isfield(C, 'n_m') || ~isfield(C, 'elevation_m') || ~isfield(C, 'x') || ~isfield(C, 'y')
        status(:) = "cross-section missing n/x/y/elevation";
        info.status = 'cross-section missing n/x/y/elevation';
        T = make_cross_section_depth_table(profileIndex, transectN_m, rowCoord, nearestRowIndex, ...
            crossSectionX_m, crossSectionY_m, crossSectionBedElevation_m, initialMedianWSE_m, crossSectionDepth_m, status);
        return
    end

    xsN = C.n_m(:);
    xsX = C.x(:);
    xsY = C.y(:);
    xsBed = C.elevation_m(:);
    goodXS = isfinite(xsN) & isfinite(xsX) & isfinite(xsY) & isfinite(xsBed);

    if nnz(goodXS) < 2
        status(:) = "fewer than two finite cross-section points";
        info.status = 'fewer than two finite cross-section points';
        T = make_cross_section_depth_table(profileIndex, transectN_m, rowCoord, nearestRowIndex, ...
            crossSectionX_m, crossSectionY_m, crossSectionBedElevation_m, initialMedianWSE_m, crossSectionDepth_m, status);
        return
    end

    xsN = xsN(goodXS);
    xsX = xsX(goodXS);
    xsY = xsY(goodXS);
    xsBed = xsBed(goodXS);

    [xsNsort, orderXS] = sort(xsN);
    xsXsort = xsX(orderXS);
    xsYsort = xsY(orderXS);
    xsBedsort = xsBed(orderXS);

    [xsNunique, ~, groupIdx] = unique(xsNsort);
    xsXunique = accumarray(groupIdx, xsXsort, [], @(v) median(v, 'omitnan'));
    xsYunique = accumarray(groupIdx, xsYsort, [], @(v) median(v, 'omitnan'));
    xsBedUnique = accumarray(groupIdx, xsBedsort, [], @(v) median(v, 'omitnan'));

    if numel(xsNunique) < 2
        status(:) = "fewer than two unique cross-section n values";
        info.status = 'fewer than two unique cross-section n values';
        T = make_cross_section_depth_table(profileIndex, transectN_m, rowCoord, nearestRowIndex, ...
            crossSectionX_m, crossSectionY_m, crossSectionBedElevation_m, initialMedianWSE_m, crossSectionDepth_m, status);
        return
    end

    crossSectionBedElevation_m = interp1(xsNunique, xsBedUnique, transectN_m(:), 'linear', NaN);
    crossSectionX_m = interp1(xsNunique, xsXunique, transectN_m(:), 'linear', NaN);
    crossSectionY_m = interp1(xsNunique, xsYunique, transectN_m(:), 'linear', NaN);

    if ~isfield(checkpoint, 'wse_map') || size(checkpoint.wse_map,1) < waterDepthWseMapAA || ...
            size(checkpoint.wse_map,2) < waterDepthWseMapIndex || ...
            isempty(checkpoint.wse_map{waterDepthWseMapAA, waterDepthWseMapIndex})
        status(:) = "water-depth WSE map unavailable";
        info.status = 'water-depth WSE map unavailable';
        T = make_cross_section_depth_table(profileIndex, transectN_m, rowCoord, nearestRowIndex, ...
            crossSectionX_m, crossSectionY_m, crossSectionBedElevation_m, initialMedianWSE_m, crossSectionDepth_m, status);
        return
    end

    Wdepth = double(checkpoint.wse_map{waterDepthWseMapAA, waterDepthWseMapIndex});
    initialWaterSurfaceElevationForDepth_m = median(Wdepth(:), 'omitnan');
    initialMedianWSE_m(:) = initialWaterSurfaceElevationForDepth_m;
    crossSectionDepth_m = initialWaterSurfaceElevationForDepth_m - crossSectionBedElevation_m;

    status(:) = "ok";
    status(~isfinite(transectN_m(:))) = "invalid transect n-location";
    status(isfinite(transectN_m(:)) & ~isfinite(crossSectionBedElevation_m)) = "transect outside cross-section n-range";
    status(isfinite(crossSectionBedElevation_m) & ~isfinite(crossSectionDepth_m)) = "depth not finite";

    info.status = 'ok';
    T = make_cross_section_depth_table(profileIndex, transectN_m, rowCoord, nearestRowIndex, ...
        crossSectionX_m, crossSectionY_m, crossSectionBedElevation_m, initialMedianWSE_m, crossSectionDepth_m, status);
end

function T = make_cross_section_depth_table(profileIndex, transectN_m, rowCoord, nearestRowIndex, ...
    crossSectionX_m, crossSectionY_m, crossSectionBedElevation_m, initialMedianWSE_m, crossSectionDepth_m, status)

    T = table( ...
        profileIndex(:), ...
        transectN_m(:), ...
        rowCoord(:), ...
        nearestRowIndex(:), ...
        crossSectionX_m(:), ...
        crossSectionY_m(:), ...
        crossSectionBedElevation_m(:), ...
        initialMedianWSE_m(:), ...
        crossSectionDepth_m(:), ...
        string(status(:)), ...
        'VariableNames', { ...
            'profileIndex', ...
            'transectN_m', ...
            'rowCoord', ...
            'nearestRowIndex', ...
            'crossSectionX_m', ...
            'crossSectionY_m', ...
            'crossSectionBedElevation_m', ...
            'initialMedianWSE_m', ...
            'crossSectionDepth_m', ...
            'crossSectionDepthStatus'});
end

function C = build_cross_section_overlay_for_depth(primaryCsvFile, fallbackCsvFiles, checkpoint, adjustToWseExtent, alignmentDirection)
% Read cross-section x/y/elevation, project to solver s/n, optionally
% straighten the cross-section along solver_n or solver_s while preserving
% original chainage and elevation order.

    C = struct('loaded', false, 'fileUsed', '');
    if nargin < 2 || isempty(fallbackCsvFiles), fallbackCsvFiles = {}; end
    if nargin < 4 || isempty(adjustToWseExtent), adjustToWseExtent = false; end
    if nargin < 5 || isempty(alignmentDirection), alignmentDirection = 'solver_n'; end

    fileUsed = resolve_existing_cross_section_file(primaryCsvFile, fallbackCsvFiles);
    if isempty(fileUsed)
        warning('Cross-section CSV not found: %s', primaryCsvFile);
        return
    end

    T = readtable(fileUsed);
    varNames = lower(strtrim(string(T.Properties.VariableNames)));
    ix = find(varNames == "x", 1, 'first');
    iy = find(varNames == "y", 1, 'first');
    iz = find(varNames == "elevation" | varNames == "z" | varNames == "wse" | ...
              varNames == "water_elevation" | varNames == "waterlevel" | ...
              varNames == "water_level", 1, 'first');

    if isempty(ix) || isempty(iy) || isempty(iz)
        error(['Cross-section CSV must contain x, y, and elevation/z columns. ', ...
               'Found columns: %s'], strjoin(T.Properties.VariableNames, ', '));
    end

    x = double(T{:, ix});
    y = double(T{:, iy});
    elevation = double(T{:, iz});
    x = x(:); y = y(:); elevation = elevation(:);

    good = isfinite(x) & isfinite(y) & isfinite(elevation);
    x = x(good); y = y(good); elevation = elevation(good);
    if isempty(x)
        warning('Cross-section CSV contains no finite x/y/elevation rows.');
        return
    end

    if ~isfield(checkpoint, 'origin') || ~isfield(checkpoint, 'e_s') || ~isfield(checkpoint, 'e_n')
        error('checkpoint.origin, checkpoint.e_s, and checkpoint.e_n are required to project cross-section x/y to solver s/n.');
    end

    origin = double(checkpoint.origin(:).');
    e_s = double(checkpoint.e_s(:));
    e_n = double(checkpoint.e_n(:));
    if numel(origin) ~= 2 || numel(e_s) ~= 2 || numel(e_n) ~= 2
        error('checkpoint.origin, checkpoint.e_s, and checkpoint.e_n must be 2-D.');
    end

    x_original = x;
    y_original = y;
    elevation_original = elevation;

    relXY = [x, y] - origin;
    s_m = relXY * e_s;
    n_m = relXY * e_n;

    adjustedToWseExtent = false;
    originalChainage_m = [0; cumsum(hypot(diff(x_original), diff(y_original)))];
    originalLength_m = originalChainage_m(end);
    adjustmentDirectionName = char(alignmentDirection);
    adjustmentDirectionSign = NaN;

    if adjustToWseExtent && isfinite(originalLength_m) && originalLength_m > 0
        switch lower(strtrim(char(alignmentDirection)))
            case {'solver_n','n','cross_channel','cross-channel','crossstream','cross-stream'}
                baseDir = e_n(:).';
                adjustmentDirectionName = 'solver_n';
            case {'solver_s','s','streamwise','along_channel','along-channel'}
                baseDir = e_s(:).';
                adjustmentDirectionName = 'solver_s';
            otherwise
                error('Unknown cross-section alignment direction: %s', alignmentDirection);
        end

        baseDir = baseDir ./ hypot(baseDir(1), baseDir(2));
        startXY = [x_original(1), y_original(1)];
        candidateSigns = [1, -1];
        candidateScores = nan(size(candidateSigns));
        candidateEndDistanceToOriginal = nan(size(candidateSigns));
        candidateX = cell(size(candidateSigns));
        candidateY = cell(size(candidateSigns));
        candidateS = cell(size(candidateSigns));
        candidateN = cell(size(candidateSigns));

        for cc = 1:numel(candidateSigns)
            thisDir = candidateSigns(cc) .* baseDir;
            xCand = startXY(1) + originalChainage_m .* thisDir(1);
            yCand = startXY(2) + originalChainage_m .* thisDir(2);
            relCand = [xCand, yCand] - origin;
            sCand = relCand * e_s;
            nCand = relCand * e_n;

            if isfield(checkpoint, 'xi') && isfield(checkpoint, 'yi')
                xiHere = checkpoint.xi(:);
                yiHere = checkpoint.yi(:);
                sTol = max([median(abs(diff(xiHere)), 'omitnan'), eps]);
                nTol = max([median(abs(diff(yiHere)), 'omitnan'), eps]);
                inside = sCand >= min(xiHere) - sTol & sCand <= max(xiHere) + sTol & ...
                         nCand >= min(yiHere) - nTol & nCand <= max(yiHere) + nTol;
                candidateScores(cc) = nnz(inside);
            else
                candidateScores(cc) = 0;
            end

            candidateEndDistanceToOriginal(cc) = hypot(xCand(end) - x_original(end), yCand(end) - y_original(end));
            candidateX{cc} = xCand;
            candidateY{cc} = yCand;
            candidateS{cc} = sCand;
            candidateN{cc} = nCand;
        end

        bestScore = max(candidateScores);
        bestIdx = find(candidateScores == bestScore);
        if numel(bestIdx) > 1
            [~, ii] = min(candidateEndDistanceToOriginal(bestIdx));
            bestIdx = bestIdx(ii);
        else
            bestIdx = bestIdx(1);
        end

        x = candidateX{bestIdx};
        y = candidateY{bestIdx};
        s_m = candidateS{bestIdx};
        n_m = candidateN{bestIdx};
        adjustmentDirectionSign = candidateSigns(bestIdx);
        adjustedToWseExtent = true;
    end

    C.loaded = true;
    C.fileUsed = fileUsed;
    C.rawTable = T;
    C.x = x;
    C.y = y;
    C.elevation_m = elevation;
    C.s_m = s_m;
    C.n_m = n_m;
    C.x_original = x_original;
    C.y_original = y_original;
    C.elevation_original_m = elevation_original;
    C.originalChainage_m = originalChainage_m;
    C.originalLength_m = originalLength_m;
    C.adjustedToWseExtent = adjustedToWseExtent;
    C.adjustmentDirectionName = adjustmentDirectionName;
    C.adjustmentDirectionSign = adjustmentDirectionSign;
    C.adjustedStartXY = [x(1), y(1)];
    C.adjustedEndXY = [x(end), y(end)];
end

function fileUsed = resolve_existing_cross_section_file(primaryCsvFile, fallbackCsvFiles)
    if nargin < 2 || isempty(fallbackCsvFiles)
        fallbackCsvFiles = {};
    end
    candidates = [{primaryCsvFile}, fallbackCsvFiles(:).'];
    fileUsed = '';
    for kk = 1:numel(candidates)
        thisFile = candidates{kk};
        if isfile(thisFile)
            fileUsed = thisFile;
            return
        end
    end
end

function functionName = choose_existing_function(candidateNames)
    functionName = '';
    for ii = 1:numel(candidateNames)
        if exist(candidateNames{ii}, 'file') == 2
            functionName = candidateNames{ii};
            return
        end
    end
    error('None of these required functions were found on the MATLAB path: %s', strjoin(candidateNames, ', '));
end

function acceptedMapIndex = resolve_accepted_map_index(checkpoint, aa, result, mapTable)
    if ~isfield(result, 'selectedMapIndex') || ~isfinite(result.selectedMapIndex)
        error('result.selectedMapIndex is missing or non-finite.');
    end

    cand = round(result.selectedMapIndex);
    if iscell(checkpoint.wse_map) && size(checkpoint.wse_map,1) >= aa && ...
            size(checkpoint.wse_map,2) >= cand && ~isempty(checkpoint.wse_map{aa,cand})
        acceptedMapIndex = cand;
        return
    end

    if ~isempty(mapTable) && istable(mapTable) && ismember('mapIndex', mapTable.Properties.VariableNames) && ...
            cand >= 1 && cand <= height(mapTable)
        acceptedMapIndex = round(mapTable.mapIndex(cand));
        return
    end

    error('Could not resolve the accepted checkpoint.wse_map index from result/mapTable.');
end

function [M, info] = get_initial_residual_map(checkpoint, aa, expectedSize)
    info = struct('source', 'none', 'index', NaN);
    M = NaN(expectedSize);

    topLevelNames = {'cell_median_val_init', 'zi', 'residualMapHist', 'residualPostMapHist', 'zi_post'};
    for ii = 1:numel(topLevelNames)
        name = topLevelNames{ii};
        if isfield(checkpoint, name)
            tmp = extract_map_from_object(checkpoint.(name), aa, 1, expectedSize);
            if ~isempty(tmp)
                M = tmp;
                info.source = ['checkpoint.', name];
                info.index = 1;
                return
            end
        end
    end

    if isfield(checkpoint, 'metrics_static')
        ms = checkpoint.metrics_static;
        metricNames = {'endIterRepResidualMapHist', 'fixedRepPathResidualMapHist', 'residualMapHist', 'residualPostMapHist'};
        for ii = 1:numel(metricNames)
            name = metricNames{ii};
            if isfield(ms, name)
                tmp = extract_map_from_object(ms.(name), aa, 1, expectedSize);
                if ~isempty(tmp)
                    M = tmp;
                    info.source = ['checkpoint.metrics_static.', name];
                    info.index = 1;
                    return
                end
            end
        end
    end

    warning('No initial residual map found. Exporting NaN initialResidualMap.');
end

function M = extract_map_from_object(obj, aa, mapIdx, expectedSize)
    M = [];
    if isempty(obj) || ~isfinite(mapIdx) || mapIdx < 1
        return
    end
    mapIdx = round(mapIdx);

    if isnumeric(obj)
        if ismatrix(obj) && isequal(size(obj), expectedSize)
            M = double(obj);
            return
        end
        if ndims(obj) == 3 && size(obj,1) == expectedSize(1) && size(obj,2) == expectedSize(2) && size(obj,3) >= mapIdx
            M = double(obj(:,:,mapIdx));
            return
        end
    end

    if iscell(obj)
        try
            if ndims(obj) >= 2 && size(obj,1) >= aa && size(obj,2) >= mapIdx && ...
                    isnumeric(obj{aa,mapIdx}) && isequal(size(obj{aa,mapIdx}), expectedSize)
                M = double(obj{aa,mapIdx});
                return
            end
        catch
        end

        try
            if numel(obj) >= mapIdx && isnumeric(obj{mapIdx}) && isequal(size(obj{mapIdx}), expectedSize)
                M = double(obj{mapIdx});
                return
            end
        catch
        end
    end
end

function S = empty_transect_struct()
    S = struct;
    S.profileIndex = NaN;
    S.rowCoord = NaN;
    S.rowCell = NaN;
    S.transectN_m = NaN;
    S.nearestRowIndex = NaN;
    S.crossSectionX_m = NaN;
    S.crossSectionY_m = NaN;
    S.crossSectionBedElevation_m = NaN;
    S.crossSectionInitialMedianWSE_m = NaN;
    S.crossSectionDepth_m = NaN;
    S.crossSectionDepthStatus = '';
    S.nIntersectedCells = NaN;
    S.nFiniteVelocityCells = NaN;
    S.initialVelocityTrackedMedian_mps = NaN;
    S.acceptedVelocityTrackedMedian_mps = NaN;
    S.observedVelocityTrackedMedian_mps = NaN;
    S.wseWavelength_m = NaN;
    S.acceptedWseAmplitude_m = NaN;
    S.acceptedWseAmplitudeStatus = '';
    S.initialResidualWavelength_m = NaN;
    S.initialResidualWavelengthStatus = '';
    S.wsePmusicWavelength_m = NaN;
    S.wsePmusicWavelengthStatus = '';
    S.wsePmusicPeakFrequency_cpm = NaN;
    S.wsePmusicPeakPower = NaN;
    S.initialResidualPmusicWavelength_m = NaN;
    S.initialResidualPmusicWavelengthStatus = '';
    S.initialResidualPmusicPeakFrequency_cpm = NaN;
    S.initialResidualPmusicPeakPower = NaN;
    S.initialResidualAmplitude = NaN;
    S.initialResidualAmplitudeUnits = '';
    S.initialResidualAmplitudeStatus = '';
    S.velocityFromWseWavelengthDeepWaterGravity_mps = NaN;
    S.velocityFromWseWavelengthDeepWaterGravityStatus = '';
    S.velocityFromResidualWavelengthDeepWaterGravity_mps = NaN;
    S.velocityFromResidualWavelengthDeepWaterGravityStatus = '';
    S.velocityFromWsePmusicWavelengthDeepWaterGravity_mps = NaN;
    S.velocityFromWsePmusicWavelengthDeepWaterGravityStatus = '';
    S.velocityFromResidualPmusicWavelengthDeepWaterGravity_mps = NaN;
    S.velocityFromResidualPmusicWavelengthDeepWaterGravityStatus = '';
    S.depthFromWseWavelength_m = NaN;
    S.depthFromWseWavelengthStatus = '';
    S.depthFromInitialResidualWavelength_m = NaN;
    S.depthFromInitialResidualWavelengthStatus = '';
    S.depthFromWsePmusicWavelength_m = NaN;
    S.depthFromWsePmusicWavelengthStatus = '';
    S.depthFromInitialResidualPmusicWavelength_m = NaN;
    S.depthFromInitialResidualPmusicWavelengthStatus = '';
    S.cellRows = [];
    S.cellCols = [];
    S.cellLinearIndices = [];
    S.velocityCellsTracked_mps = [];
end

function cells = get_transect_cells(rowCoord, colCoord, distance_m, profileIndex, mapSize)
    Ny = mapSize(1);
    Nx = mapSize(2);
    nLine = min(numel(colCoord), numel(distance_m));

    rowCell = round(rowCoord(profileIndex));
    colCells = round(colCoord(1:nLine));
    rowCells = rowCell .* ones(size(colCells));

    inGrid = isfinite(rowCells) & isfinite(colCells) & ...
        rowCells >= 1 & rowCells <= Ny & colCells >= 1 & colCells <= Nx;

    lin = sub2ind([Ny Nx], rowCells(inGrid), colCells(inGrid));
    lin = unique(lin(:), 'stable');
    [rr, cc] = ind2sub([Ny Nx], lin);

    cells = struct;
    cells.rowCell = rowCell;
    cells.rows = rr(:);
    cells.cols = cc(:);
    cells.linearIndices = lin(:);
end

function [wavelength_m, source, peakR] = get_wse_wavelength_for_profile(profileTable, acceptedMapIndex, profileIndex, acceptedWseMap, rowCoord, colCoord, distance_m, opts)
    wavelength_m = NaN;
    source = 'not found';
    peakR = NaN;

    if ~isempty(profileTable) && istable(profileTable) && ...
            all(ismember({'mapIndex','profileIndex','autocorrWavelength_m'}, profileTable.Properties.VariableNames))
        match = profileTable.mapIndex == acceptedMapIndex & profileTable.profileIndex == profileIndex;
        idx = find(match & isfinite(profileTable.autocorrWavelength_m), 1, 'first');
        if ~isempty(idx)
            wavelength_m = profileTable.autocorrWavelength_m(idx);
            source = sprintf('profileTable.autocorrWavelength_m at accepted map %d', acceptedMapIndex);
            if ismember('autocorrPeakR', profileTable.Properties.VariableNames)
                peakR = profileTable.autocorrPeakR(idx);
            end
            return
        end
    end

    R = estimate_map_wavelength_along_profile(acceptedWseMap, rowCoord, colCoord, distance_m, profileIndex, opts);
    wavelength_m = R.wavelength_m;
    source = 'acceptedWseMap sampled along transect; first positive spatial-autocorrelation peak';
    peakR = R.peakCorrelation;
end

function [amplitude_m, source, status] = get_wse_amplitude_for_profile(profileTable, acceptedMapIndex, profileIndex, acceptedWseMap, rowCoord, colCoord, distance_m, opts)
    amplitude_m = NaN;
    source = 'not found';
    status = 'not evaluated';

    if ~isempty(profileTable) && istable(profileTable) && ...
            all(ismember({'mapIndex','profileIndex','amplitude98_2_m'}, profileTable.Properties.VariableNames))
        match = profileTable.mapIndex == acceptedMapIndex & profileTable.profileIndex == profileIndex;
        idx = find(match & isfinite(profileTable.amplitude98_2_m), 1, 'first');
        if ~isempty(idx)
            amplitude_m = profileTable.amplitude98_2_m(idx);
            source = sprintf('profileTable.amplitude98_2_m at accepted map %d', acceptedMapIndex);
            status = 'ok';
            return
        end
    end

    A = estimate_map_amplitude_along_profile(acceptedWseMap, rowCoord, colCoord, distance_m, profileIndex, opts);
    amplitude_m = A.amplitude;
    source = 'acceptedWseMap sampled along transect; 0.5*(P98-P2) after linear detrending';
    status = A.status;
end

function A = estimate_map_amplitude_along_profile(M, rowCoord, colCoord, distance_m, profileIndex, opts)
    A = struct;
    A.amplitude = NaN;
    A.units = 'map units';
    A.status = 'not evaluated';
    A.profileIndex = profileIndex;
    A.rowCoord = rowCoord(profileIndex);
    A.segmentIndices = [];
    A.segmentDistance_m = [];
    A.segmentValues = [];
    A.segmentDetrendedValues = [];
    A.amplitudeSource = 'map sampled along transect; 0.5*(P98-P2) after linear detrending';

    if isfield(opts, 'ampLowPct') && isfinite(opts.ampLowPct)
        pLo = opts.ampLowPct;
    else
        pLo = 2;
    end
    if isfield(opts, 'ampHighPct') && isfinite(opts.ampHighPct)
        pHi = opts.ampHighPct;
    else
        pHi = 98;
    end

    nLine = min(numel(colCoord), numel(distance_m));
    rowLine = rowCoord(profileIndex) .* ones(nLine, 1);
    colLine = colCoord(1:nLine);
    x = distance_m(1:nLine);
    z = interp2(M, colLine, rowLine, 'linear', NaN);

    finiteIdx = find(isfinite(x) & isfinite(z));
    if numel(finiteIdx) < 5
        A.status = 'too few finite samples';
        return
    end

    seg = finiteIdx(1):finiteIdx(end);
    xSeg = x(seg);
    zSeg = z(seg);
    xSeg = xSeg - xSeg(1);

    zDet = detrend_profile_linear_local(xSeg, zSeg);
    A.amplitude = robust_half_range_local(zDet, pLo, pHi);
    A.segmentIndices = seg(:);
    A.segmentDistance_m = xSeg(:);
    A.segmentValues = zSeg(:);
    A.segmentDetrendedValues = zDet(:);

    if isfinite(A.amplitude)
        A.status = 'ok';
    else
        A.status = 'amplitude not finite';
    end
end

function zDet = detrend_profile_linear_local(x, z)
    x = x(:);
    z = z(:);
    zDet = NaN(size(z));
    good = isfinite(x) & isfinite(z);
    if nnz(good) < 3
        return
    end
    xv = x(good);
    zv = z(good);
    xs = (xv - mean(xv)) ./ max(std(xv), eps);
    p = polyfit(xs, zv, 1);
    xsAll = (x - mean(xv)) ./ max(std(xv), eps);
    trend = polyval(p, xsAll);
    zDet = z - trend;
    zDet(good) = zDet(good) - mean(zDet(good));
    zDet(~good) = NaN;
end

function a = robust_half_range_local(v, pLo, pHi)
    v = v(isfinite(v));
    if numel(v) < 5
        a = NaN;
    else
        a = 0.5 .* (prctile(v, pHi) - prctile(v, pLo));
    end
end

function R = estimate_map_wavelength_along_profile(M, rowCoord, colCoord, distance_m, profileIndex, opts)
    nLine = min(numel(colCoord), numel(distance_m));
    rowLine = rowCoord(profileIndex) .* ones(nLine, 1);
    colLine = colCoord(1:nLine);
    x = distance_m(1:nLine);

    z = interp2(M, colLine, rowLine, 'linear', NaN);

    finiteIdx = find(isfinite(x) & isfinite(z));
    if numel(finiteIdx) >= 2
        seg = finiteIdx(1):finiteIdx(end);
    else
        seg = 1:nLine;
    end

    xSeg = x(seg);
    zSeg = z(seg);
    xSeg = xSeg - xSeg(1);

    R = estimate_wavelength_from_spatial_autocorrelation_local(xSeg, zSeg, opts);
    R.profileIndex = profileIndex;
    R.rowCoord = rowCoord(profileIndex);
    R.segmentIndices = seg(:);
    R.segmentDistance_m = xSeg(:);
    R.segmentValues = zSeg(:);
    R.wavelengthSource = 'map sampled along transect; first positive spatial-autocorrelation peak';
end

function R = estimate_map_wavelength_pmusic_along_profile(M, rowCoord, colCoord, distance_m, profileIndex, opts)
    nLine = min(numel(colCoord), numel(distance_m));
    rowLine = rowCoord(profileIndex) .* ones(nLine, 1);
    colLine = colCoord(1:nLine);
    x = distance_m(1:nLine);

    z = interp2(M, colLine, rowLine, 'linear', NaN);

    finiteIdx = find(isfinite(x) & isfinite(z));
    if numel(finiteIdx) >= 2
        seg = finiteIdx(1):finiteIdx(end);
    else
        seg = 1:nLine;
    end

    xSeg = x(seg);
    zSeg = z(seg);
    xSeg = xSeg - xSeg(1);

    R = estimate_wavelength_from_pmusic_local(xSeg, zSeg, opts);
    R.profileIndex = profileIndex;
    R.rowCoord = rowCoord(profileIndex);
    R.segmentIndices = seg(:);
    R.segmentDistance_m = xSeg(:);
    R.segmentValues = zSeg(:);
    R.wavelengthSource = 'map sampled along transect; PMUSIC pseudospectrum peak';
end

function R = estimate_wavelength_from_pmusic_local(x, z, opts)
    if nargin < 3 || isempty(opts), opts = struct; end
    opts = default_struct_field(opts, 'usePmusicWavelength', true);
    opts = default_struct_field(opts, 'pmusicModelOrder', 2);
    opts = default_struct_field(opts, 'pmusicNfft', 2048);
    opts = default_struct_field(opts, 'pmusicWindowLength', []);
    opts = default_struct_field(opts, 'pmusicOverlap', []);
    opts = default_struct_field(opts, 'minWavelength_m', 0.30);
    opts = default_struct_field(opts, 'maxWavelength_m', Inf);
    opts = default_struct_field(opts, 'maxLagFraction', 0.75);
    opts = default_struct_field(opts, 'detrendOrder', 1);

    R = struct('wavelength_m',NaN, ...
               'peakFrequency_cpm',NaN, ...
               'peakPower',NaN, ...
               'status','not evaluated');

    if ~opts.usePmusicWavelength
        R.status = 'pmusic disabled by opts.usePmusicWavelength';
        return
    end
    if exist('pmusic', 'file') ~= 2
        R.status = 'pmusic unavailable: Signal Processing Toolbox function not found';
        return
    end

    x = x(:);
    z = z(:);
    valid = isfinite(x) & isfinite(z);
    x = x(valid);
    z = z(valid);

    if numel(x) < 12
        R.status = 'too few samples';
        return
    end

    [x, order] = sort(x);
    z = z(order);
    [x, uniqueIdx] = unique(x, 'stable');
    z = z(uniqueIdx);

    dx = median(diff(x), 'omitnan');
    if ~isfinite(dx) || dx <= 0
        R.status = 'invalid dx';
        return
    end

    xUniform = (min(x):dx:max(x)).';
    zUniform = interp1(x, z, xUniform, 'linear', 'extrap');
    n = numel(zUniform);

    if n < 12
        R.status = 'too short';
        return
    end

    if opts.detrendOrder >= 0 && n > opts.detrendOrder + 1
        pDet = polyfit(xUniform, zUniform, opts.detrendOrder);
        zProcessed = zUniform - polyval(pDet, xUniform);
    else
        zProcessed = zUniform;
    end

    zProcessed = zProcessed - mean(zProcessed, 'omitnan');
    if sum(zProcessed.^2, 'omitnan') <= eps
        R.status = 'no variation';
        return
    end

    fsSpatial = 1 ./ dx; % samples per metre, so PMUSIC frequencies are cycles per metre
    modelOrder = max(1, round(opts.pmusicModelOrder));
    nfft = max(16, round(opts.pmusicNfft));

    if isempty(opts.pmusicWindowLength) || ~isfinite(opts.pmusicWindowLength)
        nwin = min(n, max(modelOrder + 2, floor(n/2)));
    else
        nwin = min(n, max(modelOrder + 2, round(opts.pmusicWindowLength)));
    end

    if isempty(opts.pmusicOverlap) || ~isfinite(opts.pmusicOverlap)
        noverlap = floor(0.5 .* nwin);
    else
        noverlap = max(0, min(nwin - 1, round(opts.pmusicOverlap)));
    end

    if nwin <= modelOrder
        R.status = 'pmusic window length must be greater than model order';
        return
    end

    try
        [Pmusic, f] = pmusic(zProcessed, modelOrder, nfft, fsSpatial, nwin, noverlap);
    catch ME
        R.status = ['pmusic failed: ', ME.message];
        return
    end

    Pmusic = Pmusic(:);
    f = f(:);

    profileLength_m = max(xUniform) - min(xUniform);
    maxWavelengthAllowed = opts.maxWavelength_m;
    if ~isfinite(maxWavelengthAllowed)
        maxWavelengthAllowed = max(opts.minWavelength_m, opts.maxLagFraction .* profileLength_m);
    end

    fMin = 1 ./ maxWavelengthAllowed;
    fMax = 1 ./ opts.minWavelength_m;
    good = isfinite(Pmusic) & isfinite(f) & f > 0 & f >= fMin & f <= fMax;

    if ~any(good)
        R.status = 'no PMUSIC frequencies inside wavelength bounds';
        return
    end

    idxGood = find(good);
    [peakPower, localIdx] = max(Pmusic(good));
    idxPeak = idxGood(localIdx);
    fPeak = f(idxPeak);

    if ~isfinite(fPeak) || fPeak <= 0
        R.status = 'invalid PMUSIC peak frequency';
        return
    end

    R.wavelength_m = 1 ./ fPeak;
    R.peakFrequency_cpm = fPeak;
    R.peakPower = peakPower;
    R.status = 'ok';
end

function V = velocity_from_deep_water_gravity_wave(lambda_m)
    V = struct;
    V.velocity_mps = NaN;
    V.wavelength_m = lambda_m;
    V.status = 'not evaluated';

    if ~isfinite(lambda_m) || lambda_m <= 0
        V.status = 'invalid wavelength';
        return
    end

    g = 9.81;
    V.velocity_mps = sqrt(g .* lambda_m ./ (2*pi));
    V.status = 'ok';
end


function depthEstimate = estimate_depth_from_known_wavelength_filtered(vel, knownWavelength_m, opts)
    % Estimate depth from observed velocity and known wavelength.
    % This uses the observed dt=0.5 s tracked velocity plus a measured
    % wavelength source to estimate depth.
    depthEstimate = struct;
    depthEstimate.velocityInput_mps = vel;
    depthEstimate.knownWavelength_m = knownWavelength_m;
    depthEstimate.depthEstimate_m = NaN;
    depthEstimate.rawDepthEstimate_m = NaN;
    depthEstimate.wavelengthAtEstimate_m = NaN;
    depthEstimate.wavelengthMismatchFraction = NaN;
    depthEstimate.nearestIndex = NaN;
    depthEstimate.depthGrid_m = opts.depthGrid_m;
    depthEstimate.wavelengthCurve_m = NaN(size(opts.depthGrid_m));
    depthEstimate.k0 = NaN(size(opts.depthGrid_m));
    depthEstimate.status = 'not evaluated';

    if ~isfinite(vel) || vel <= 0
        depthEstimate.status = 'discarded: invalid 0.5 s velocity input';
        return
    end
    if ~isfinite(knownWavelength_m) || knownWavelength_m <= 0
        depthEstimate.status = 'discarded: invalid wavelength input';
        return
    end

    g = 9.81;
    surftens = 72.75e-03; % N/m, approx. 20 deg C
    dens = 998.2;         % kg/m^3, approx. 20 deg C
    alpha = 0.83;
    m = 2*(1-alpha);

    depth = opts.depthGrid_m;
    k0init = 2*pi;

    fun = @(k,d) (m + ((g*d)/vel^2) * ...
        (1 + dens*g/(surftens*k^2)) / (dens*g/(surftens*k^2))) * ...
        tanh(k*d)/(k*d) - 1;

    k0 = NaN(size(depth));
    for idepth = 1:length(depth)
        d = depth(idepth);
        if vel^2./(g*d) <= 1
            try
                k0(idepth) = fzero(@(k) fun(k,d), k0init);
            catch
                k0(idepth) = NaN;
            end
        else
            k0(idepth) = NaN;
        end
    end

    wavelength = 2*pi ./ k0;

    depthEstimate.depthGrid_m = depth;
    depthEstimate.wavelengthCurve_m = wavelength;
    depthEstimate.k0 = k0;

    finiteCurve = isfinite(wavelength) & isfinite(depth);
    if ~any(finiteCurve)
        depthEstimate.status = 'discarded: no finite wavelength-depth curve values';
        return
    end

    finiteWavelengths = wavelength(finiteCurve);
    if knownWavelength_m < min(finiteWavelengths) || knownWavelength_m > max(finiteWavelengths)
        depthEstimate.status = 'discarded: wavelength outside finite curve range';
        return
    end

    idxFinite = find(finiteCurve);
    [~, localIdx] = min(abs(wavelength(finiteCurve) - knownWavelength_m));
    nearestIdx = idxFinite(localIdx);

    rawDepth = depth(nearestIdx);
    wavelengthAtEstimate = wavelength(nearestIdx);
    mismatchFraction = abs(wavelengthAtEstimate - knownWavelength_m) ./ max(abs(knownWavelength_m), eps);

    depthEstimate.rawDepthEstimate_m = rawDepth;
    depthEstimate.wavelengthAtEstimate_m = wavelengthAtEstimate;
    depthEstimate.wavelengthMismatchFraction = mismatchFraction;
    depthEstimate.nearestIndex = nearestIdx;

    if mismatchFraction > opts.maxWavelengthMismatchFraction
        depthEstimate.status = 'discarded: nearest curve wavelength mismatch too large';
        return
    end
    if rawDepth < opts.reasonableDepthRange_m(1) || rawDepth > opts.reasonableDepthRange_m(2)
        depthEstimate.status = 'discarded: depth outside reasonable range';
        return
    end

    depthEstimate.depthEstimate_m = rawDepth;
    depthEstimate.status = 'ok';
end

function R = estimate_wavelength_from_spatial_autocorrelation_local(x, z, opts)
    if nargin < 3 || isempty(opts), opts = struct; end
    opts = default_struct_field(opts, 'minWavelength_m', 0.30);
    opts = default_struct_field(opts, 'maxWavelength_m', Inf);
    opts = default_struct_field(opts, 'maxLagFraction', 0.75);
    opts = default_struct_field(opts, 'minPeakCorrelation', 0.10);
    opts = default_struct_field(opts, 'detrendOrder', 1);

    R = struct('wavelength_m',NaN,'peakCorrelation',NaN,'peakLag_m',NaN,'status','not evaluated');

    x = x(:);
    z = z(:);
    valid = isfinite(x) & isfinite(z);
    x = x(valid);
    z = z(valid);

    if numel(x) < 8
        R.status = 'too few samples';
        return
    end

    [x, order] = sort(x);
    z = z(order);
    [x, uniqueIdx] = unique(x, 'stable');
    z = z(uniqueIdx);

    dx = median(diff(x), 'omitnan');
    if ~isfinite(dx) || dx <= 0
        R.status = 'invalid dx';
        return
    end

    xUniform = (min(x):dx:max(x)).';
    zUniform = interp1(x, z, xUniform, 'linear', 'extrap');
    n = numel(zUniform);

    if n < 8
        R.status = 'too short';
        return
    end

    if opts.detrendOrder >= 0 && n > opts.detrendOrder + 1
        p = polyfit(xUniform, zUniform, opts.detrendOrder);
        zProcessed = zUniform - polyval(p, xUniform);
    else
        zProcessed = zUniform;
    end

    zProcessed = zProcessed - mean(zProcessed, 'omitnan');
    if sum(zProcessed.^2, 'omitnan') <= eps
        R.status = 'no variation';
        return
    end

    maxLagByFraction = floor(max(0.05, min(opts.maxLagFraction, 0.95)) * (n - 1));
    maxLagByWavelength = n - 2;
    if isfinite(opts.maxWavelength_m)
        maxLagByWavelength = floor(opts.maxWavelength_m / dx);
    end
    maxLag = min([n-2, maxLagByFraction, maxLagByWavelength]);
    minLag = max(1, ceil(opts.minWavelength_m / dx));

    if maxLag < minLag + 2
        R.status = 'lag range too short';
        return
    end

    acf = NaN(maxLag+1,1);
    for lag = 0:maxLag
        a = zProcessed(1:n-lag);
        b = zProcessed(1+lag:n);
        den = sqrt(sum(a.^2, 'omitnan') * sum(b.^2, 'omitnan'));
        if den > eps
            acf(lag+1) = sum(a.*b, 'omitnan') / den;
        end
    end

    idxPeak = [];
    for ii = minLag+1:maxLag
        if ii > 1 && ii < numel(acf) && isfinite(acf(ii)) && ...
                acf(ii) >= opts.minPeakCorrelation && acf(ii) > acf(ii-1) && acf(ii) >= acf(ii+1)
            idxPeak = ii;
            break
        end
    end

    if isempty(idxPeak)
        R.status = 'no clear peak';
        return
    end

    delta = 0;
    if idxPeak > 1 && idxPeak < numel(acf) && all(isfinite(acf([idxPeak-1 idxPeak idxPeak+1])))
        den = acf(idxPeak-1) - 2*acf(idxPeak) + acf(idxPeak+1);
        if abs(den) > eps
            delta = max(-1, min(1, 0.5*(acf(idxPeak-1)-acf(idxPeak+1))/den));
        end
    end

    lagSamples = (idxPeak-1) + delta;
    R.wavelength_m = lagSamples * dx;
    R.peakCorrelation = acf(idxPeak);
    R.peakLag_m = R.wavelength_m;
    R.status = 'ok';
end

function s = default_struct_field(s, f, v)
    if ~isfield(s, f) || isempty(s.(f))
        s.(f) = v;
    end
end

function Tstruct = remove_large_fields_for_table(S)
    Tstruct = rmfield(S, {'cellRows','cellCols','cellLinearIndices','velocityCellsTracked_mps'});
end
