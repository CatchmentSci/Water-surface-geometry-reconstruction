function [result, mapTable, profileTable] = KLT_select_wave_checkpoint_refined(checkpoint, opts)
%KLT_SELECT_WAVE_CHECKPOINT_REFINED_FIRST10M
% Best-available checkpoint-only stop and wave-parameter estimator.
%
% This function uses only data already saved inside the KLT solver checkpoint:
%   - checkpoint.wse_map
%   - checkpoint.metrics_static.fixedRepPathMedianAbsResidualHist
%   - checkpoint.metrics_static.endIterRepCellMedianAbsResidualHist
%   - checkpoint.metrics_static.localMedianRelativeImprovementHist
%   - checkpoint.metrics_static surface/update diagnostics when available
%
% It does not require known amplitude/wavelength and does not assume a sine wave.
% Wavelength is estimated from the first positive spatial-autocorrelation peak
% and, additionally, from a PMUSIC spatial pseudospectrum peak.
% Amplitude is estimated as 0.5*(P98-P2) of detrended transect profiles.
%
% Example:
%   S = load('your_checkpoint.mat','checkpoint');
%   [result, mapTable, profileTable] = KLT_select_wave_checkpoint_refined_first10m(S.checkpoint);
%
% Result fields:
%   selectedMapIndex, stopMethod, amplitude_m, wavelength_m,
%   amplitudeIQR_m, wavelengthIQR_m, nAmplitudeProfiles, nWavelengthProfiles,
%   pmusicWavelength_m, pmusicWavelengthIQR_m, nPmusicWavelengthProfiles

    if nargin < 2 || isempty(opts)
        opts = struct;
    end

    % Compatibility: allow callers to pass the plotter-style opts struct
    % directly, e.g. opts.selectorOpts.cameraRowOrder = 'ascending'.
    % Fields inside selectorOpts are copied into opts unless opts already
    % defines the same selector option at the top level.
    if isstruct(opts) && isfield(opts, 'selectorOpts') && isstruct(opts.selectorOpts)
        nestedSelectorOpts = opts.selectorOpts;
        nestedFields = fieldnames(nestedSelectorOpts);
        for ff = 1:numel(nestedFields)
            thisField = nestedFields{ff};
            if ~isfield(opts, thisField) || isempty(opts.(thisField))
                opts.(thisField) = nestedSelectorOpts.(thisField);
            end
        end
    end

    opts = set_default_opts(opts);

    if ~isfield(checkpoint, 'wse_map')
        error('checkpoint.wse_map is required.');
    end
    if ~isfield(checkpoint, 'xi') || ~isfield(checkpoint, 'yi')
        error('checkpoint.xi and checkpoint.yi are required.');
    end

    if isfield(checkpoint, 'aa') && isfinite(checkpoint.aa)
        aa = checkpoint.aa;
    else
        aa = 1;
    end

    [wseMaps, mapIndices] = collect_wse_maps(checkpoint.wse_map, aa);
    nMaps = numel(wseMaps);
    if nMaps < 4
        error('At least four non-empty WSE maps are required.');
    end

    xi = checkpoint.xi(:).';
    yi = checkpoint.yi(:);
    [Ny, Nx] = size(wseMaps{1});

    if numel(xi) ~= Nx || numel(yi) ~= Ny
        error('checkpoint.xi/yi sizes do not match WSE map dimensions.');
    end

    referenceMask = build_reference_mask(wseMaps{1}, wseMaps{end}, opts.backgroundTol, opts.minBlobSizeCells);
    [rowCoord, colCoord, distance_m, stabilityBand] = build_profile_geometry(xi, yi, referenceMask, opts, checkpoint);
    nProfiles = numel(rowCoord);

    A = nan(nMaps, nProfiles);
    L = nan(nMaps, nProfiles);
    R = nan(nMaps, nProfiles);
    Lpmusic = nan(nMaps, nProfiles);
    Fpmusic = nan(nMaps, nProfiles);
    Ppmusic = nan(nMaps, nProfiles);
    pmusicStatus = repmat("not evaluated", nMaps, nProfiles);
    nFinite = nan(nMaps, nProfiles);
    profileShape = cell(nMaps, nProfiles);
    relUpdate = nan(nMaps, 1);
    mapRelief = nan(nMaps, 1);

    prevW = [];
    for mm = 1:nMaps
        W = wseMaps{mm};
        W(~referenceMask) = NaN;

        mapRelief(mm) = robust_half_range(W(referenceMask), opts.ampLowPct, opts.ampHighPct);
        if ~isempty(prevW)
            dW = W - prevW;
            dA = robust_half_range(dW(referenceMask), opts.ampLowPct, opts.ampHighPct);
            relUpdate(mm) = dA ./ max(mapRelief(mm), eps);
        end
        prevW = W;

        for pp = 1:nProfiles
            z = sample_profile(W, rowCoord(pp), colCoord);
            finiteIdx = find(isfinite(z));
            if numel(finiteIdx) < opts.minFiniteSamples
                continue
            end

            seg = finiteIdx(1):finiteIdx(end);
            zSeg = nan(size(z));
            zSeg(seg) = z(seg);

            zDet = robust_linear_detrend(distance_m, zSeg);
            A(mm, pp) = robust_half_range(zDet, opts.ampLowPct, opts.ampHighPct);
            nFinite(mm, pp) = nnz(isfinite(zDet));

            if isfinite(A(mm, pp)) && A(mm, pp) > opts.amplitudeNoiseFloor_m
                profileShape{mm, pp} = winsorise(zDet ./ A(mm, pp), 1, 99);
            else
                profileShape{mm, pp} = nan(size(zDet));
            end

            maxLag_m = min(opts.maxWavelength_m, 0.75 * (distance_m(seg(end)) - distance_m(seg(1))));
            [L(mm, pp), R(mm, pp)] = autocorr_wavelength(distance_m, zSeg, opts.minWavelength_m, maxLag_m, opts.minAutocorrPeakR);

            Pm = pmusic_wavelength(distance_m, zSeg, opts.minWavelength_m, maxLag_m, opts);
            Lpmusic(mm, pp) = Pm.wavelength_m;
            Fpmusic(mm, pp) = Pm.peakFrequency_cpm;
            Ppmusic(mm, pp) = Pm.peakPower;
            pmusicStatus(mm, pp) = string(Pm.status);
        end
    end

    [Amed, AiqrLow, AiqrHigh, nA] = aggregate_profiles(A, true(size(A)), opts.minProfiles);
    [Lmed, LiqrLow, LiqrHigh, nL] = aggregate_profiles(L, R >= opts.minAutocorrPeakR, opts.minProfiles);
    [LpmMed, LpmIqrLow, LpmIqrHigh, nLpm] = aggregate_profiles(Lpmusic, isfinite(Lpmusic), opts.minProfiles);

    ampIqrFrac = (AiqrHigh - AiqrLow) ./ max(Amed, eps);
    wavIqrFrac = (LiqrHigh - LiqrLow) ./ max(Lmed, eps);

    shapeChange = nan(nMaps, 1);
    for mm = 2:nMaps
        d = nan(nProfiles, 1);
        for pp = 1:nProfiles
            a = profileShape{mm, pp};
            b = profileShape{mm-1, pp};
            good = isfinite(a) & isfinite(b);
            if nnz(good) >= opts.minFiniteSamples
                d(pp) = robust_rms(a(good) - b(good));
            end
        end
        goodD = isfinite(d);
        if nnz(goodD) >= opts.minProfiles
            shapeChange(mm) = median(d(goodD), 'omitnan');
        end
    end

    M = metrics_to_map_series(checkpoint, nMaps);
    [selectedMap, stopMethod, stopInfo] = select_stop_map(Amed, Lmed, nA, nL, ampIqrFrac, wavIqrFrac, shapeChange, relUpdate, M, opts);

    mapTable = table( ...
        mapIndices(:), Amed(:), AiqrLow(:), AiqrHigh(:), nA(:), ...
        Lmed(:), LiqrLow(:), LiqrHigh(:), nL(:), ...
        LpmMed(:), LpmIqrLow(:), LpmIqrHigh(:), nLpm(:), ...
        ampIqrFrac(:), wavIqrFrac(:), shapeChange(:), relUpdate(:), mapRelief(:), ...
        M.localMedianRelativeImprovement(:), M.fixedRepMedianAbsResidual(:), M.endRepCellMedianAbsResidual(:), ...
        'VariableNames', { ...
            'mapIndex','amplitudeMedian_m','amplitudeQ25_m','amplitudeQ75_m','nAmplitudeProfiles', ...
            'wavelengthMedian_m','wavelengthQ25_m','wavelengthQ75_m','nWavelengthProfiles', ...
            'pmusicWavelengthMedian_m','pmusicWavelengthQ25_m','pmusicWavelengthQ75_m','nPmusicWavelengthProfiles', ...
            'amplitudeIQRFraction','wavelengthIQRFraction','shapeChange','relativeMapUpdate','mapRelief_m', ...
            'localMedianRelativeImprovement','fixedRepMedianAbsResidual','endRepCellMedianAbsResidual'});

    profileRows = [];
    profileCols = [];
    profileA = [];
    profileL = [];
    profileR = [];
    profileLpm = [];
    profileFpm = [];
    profilePpm = [];
    profilePmStatus = strings(0, 1);
    for mm = 1:nMaps
        profileRows = [profileRows; repmat(mapIndices(mm), nProfiles, 1)]; %#ok<AGROW>
        profileCols = [profileCols; (1:nProfiles).']; %#ok<AGROW>
        profileA = [profileA; A(mm,:).']; %#ok<AGROW>
        profileL = [profileL; L(mm,:).']; %#ok<AGROW>
        profileR = [profileR; R(mm,:).']; %#ok<AGROW>
        profileLpm = [profileLpm; Lpmusic(mm,:).']; %#ok<AGROW>
        profileFpm = [profileFpm; Fpmusic(mm,:).']; %#ok<AGROW>
        profilePpm = [profilePpm; Ppmusic(mm,:).']; %#ok<AGROW>
        profilePmStatus = [profilePmStatus; pmusicStatus(mm,:).']; %#ok<AGROW>
    end
    profileTable = table(profileRows, profileCols, profileA, profileL, profileR, profileLpm, profileFpm, profilePpm, profilePmStatus, ...
        'VariableNames', {'mapIndex','profileIndex','amplitude98_2_m','autocorrWavelength_m','autocorrPeakR', ...
                          'pmusicWavelength_m','pmusicPeakFrequency_cpm','pmusicPeakPower','pmusicStatus'});

    selectedMapRow = selectedMap;
    selectedMapIndex = mapIndices(selectedMapRow);

    result = struct;
    result.selectedMapRow = selectedMapRow;
    result.selectedMapIndex = selectedMapIndex;
    result.stopMethod = stopMethod;
    result.stopInfo = stopInfo;
    result.amplitude_m = mapTable.amplitudeMedian_m(selectedMapRow);
    result.amplitudeIQR_m = [mapTable.amplitudeQ25_m(selectedMapRow), mapTable.amplitudeQ75_m(selectedMapRow)];
    result.wavelength_m = mapTable.wavelengthMedian_m(selectedMapRow);
    result.wavelengthIQR_m = [mapTable.wavelengthQ25_m(selectedMapRow), mapTable.wavelengthQ75_m(selectedMapRow)];
    result.pmusicWavelength_m = mapTable.pmusicWavelengthMedian_m(selectedMapRow);
    result.pmusicWavelengthIQR_m = [mapTable.pmusicWavelengthQ25_m(selectedMapRow), mapTable.pmusicWavelengthQ75_m(selectedMapRow)];
    result.nAmplitudeProfiles = mapTable.nAmplitudeProfiles(selectedMapRow);
    result.nWavelengthProfiles = mapTable.nWavelengthProfiles(selectedMapRow);
    result.nPmusicWavelengthProfiles = mapTable.nPmusicWavelengthProfiles(selectedMapRow);
    result.referenceMask = referenceMask;
    result.rowCoord = rowCoord;
    result.colCoord = colCoord;
    result.distance_m = distance_m;
    result.stabilityBand = stabilityBand;
    result.stabilityAssessmentMode = opts.stabilityAssessmentMode;
    result.firstPopulatedRowIndex = stabilityBand.firstPopulatedRowIndex;
    result.firstPopulatedRowCoord_m = stabilityBand.firstPopulatedRowCoord_m;
    result.stabilityBandLength_m = stabilityBand.length_m;
    result.stabilityBandRows = stabilityBand.rows(:);
    result.waveParameterDisplay = table( ...
        result.selectedMapIndex, result.amplitude_m, result.wavelength_m, result.pmusicWavelength_m, ...
        result.nWavelengthProfiles, result.nPmusicWavelengthProfiles, ...
        'VariableNames', {'selectedMapIndex','amplitude_m','autocorrWavelength_m','pmusicWavelength_m', ...
                          'nAutocorrWavelengthProfiles','nPmusicWavelengthProfiles'});

    if opts.displaySummary
        fprintf('\nKLT_select_wave_checkpoint_refined selected map %d (%s)\n', selectedMapIndex, stopMethod);
        fprintf('  Amplitude median: %.5g m, n = %d profiles\n', result.amplitude_m, result.nAmplitudeProfiles);
        fprintf('  Autocorr wavelength median: %.5g m, IQR [%.5g %.5g] m, n = %d profiles\n', ...
            result.wavelength_m, result.wavelengthIQR_m(1), result.wavelengthIQR_m(2), result.nWavelengthProfiles);
        fprintf('  PMUSIC wavelength median: %.5g m, IQR [%.5g %.5g] m, n = %d profiles\n', ...
            result.pmusicWavelength_m, result.pmusicWavelengthIQR_m(1), result.pmusicWavelengthIQR_m(2), result.nPmusicWavelengthProfiles);
        disp(result.waveParameterDisplay);
    end
end

function opts = set_default_opts(opts)
    % Stability-assessment profile window.
    % Default is now camera-relative: find the first WSE-populated row
    % from the camera side and analyse only the next firstSensedLength_m
    % moving away from the camera.
    %
    % Set opts.stabilityAssessmentMode = 'manual_row_range' to recover the
    % older fixed-row behaviour using opts.rowRange.
    opts = default_field(opts, 'stabilityAssessmentMode', 'first_populated_10m_from_camera');
    opts = default_field(opts, 'firstSensedLength_m', 10);
    opts = default_field(opts, 'cameraRowOrder', 'auto'); % 'auto', 'descending', or 'ascending'
    opts = default_field(opts, 'minPopulatedCellsPerRow', []);
    opts = default_field(opts, 'rowRange', [180 360]);
    opts = default_field(opts, 'profileSpacing_m', 0.5);
    opts = default_field(opts, 'sampleSpacingAlongProfile_m', []);
    opts = default_field(opts, 'backgroundTol', 1e-6);
    opts = default_field(opts, 'minBlobSizeCells', 20);
    opts = default_field(opts, 'minFiniteSamples', 12);
    opts = default_field(opts, 'minProfiles', 8);
    opts = default_field(opts, 'ampLowPct', 2);
    opts = default_field(opts, 'ampHighPct', 98);
    opts = default_field(opts, 'amplitudeNoiseFloor_m', 0.003);
    opts = default_field(opts, 'minWavelength_m', 0.30);
    opts = default_field(opts, 'maxWavelength_m', Inf);
    opts = default_field(opts, 'minAutocorrPeakR', 0.10);
    opts = default_field(opts, 'usePmusicWavelength', true);
    opts = default_field(opts, 'pmusicModelOrder', 2);          % 2 for one real sinusoidal component
    opts = default_field(opts, 'pmusicNfft', 2048);
    opts = default_field(opts, 'pmusicWindowLength', []);       % [] chooses a data-dependent value
    opts = default_field(opts, 'pmusicOverlap', []);            % [] chooses 50% overlap
    opts = default_field(opts, 'displaySummary', true);
    opts = default_field(opts, 'weakUpdateRelMax', 0.12);
    opts = default_field(opts, 'weakShapeTol', 0.14);
    opts = default_field(opts, 'weakAmpIqrFracMax', 0.45);
    opts = default_field(opts, 'weakWavIqrFracMax', 0.25);
    opts = default_field(opts, 'wavStabRel', 0.05);
    opts = default_field(opts, 'wavStabAbs_m', 0.05);
    opts = default_field(opts, 'stableWindow', 3);
    opts = default_field(opts, 'localImprovementCollapseThreshold', 0.05);
    opts = default_field(opts, 'residualNearMinimumFraction', 0.20);
    opts = default_field(opts, 'strongForwardWindow', 5);
    opts = default_field(opts, 'strongMaxForwardAmpGrowth', 0.05);
    opts = default_field(opts, 'strongMaxForwardAmpAbsGrowth_m', 0.003);
    opts = default_field(opts, 'weakForwardWindow', 5);
    opts = default_field(opts, 'weakMaxForwardAmpGrowth', 0.05);
    opts = default_field(opts, 'weakMaxForwardAmpAbsGrowth_m', 0.003);
    opts = default_field(opts, 'requireForwardAmplitudePlateau', true);
    opts = default_field(opts, 'noPlateauFallbackMode', 'latest_supported');
    opts = default_field(opts, 'strongShapeTol', 0.14);
    opts = default_field(opts, 'strongAmpIqrFracMax', 0.45);
    opts = default_field(opts, 'strongWavIqrFracMax', 0.25);
end

function s = default_field(s, f, v)
    if ~isfield(s, f) || isempty(s.(f))
        s.(f) = v;
    end
end

function [maps, mapIndices] = collect_wse_maps(wse_map, aa)
    maps = {};
    mapIndices = [];
    if iscell(wse_map)
        for kk = 1:size(wse_map, 2)
            W = wse_map{aa, kk};
            if ~isempty(W) && isnumeric(W) && ismatrix(W)
                maps{end+1,1} = double(W); %#ok<AGROW>
                mapIndices(end+1,1) = kk; %#ok<AGROW>
            end
        end
    else
        error('checkpoint.wse_map must be a cell array.');
    end
end

function mask = build_reference_mask(W0, Wref, tol, minBlob)
    base = median(W0(isfinite(W0)), 'omitnan');
    mask0 = isfinite(Wref) & abs(Wref - base) > tol;
    try
        CC = bwconncomp(mask0);
        if CC.NumObjects > 0
            sizes = cellfun(@numel, CC.PixelIdxList);
            keep = find(sizes >= minBlob);
            if isempty(keep)
                [~, k] = max(sizes);
            else
                [~, kk] = max(sizes(keep));
                k = keep(kk);
            end
            mask = false(size(mask0));
            mask(CC.PixelIdxList{k}) = true;
            mask = imfill(mask, 'holes');
        else
            mask = mask0;
        end
    catch
        mask = mask0;
    end
end

function [rowCoord, colCoord, distance_m, stabilityBand] = build_profile_geometry(xi, yi, mask, opts, checkpoint)
    Ny = numel(yi);
    Nx = numel(xi);

    mode = lower(char(opts.stabilityAssessmentMode));

    if strcmp(mode, 'manual_row_range') || strcmp(mode, 'rowrange') || strcmp(mode, 'manual')
        rowA = max(1, min(Ny, round(opts.rowRange(1))));
        rowB = max(1, min(Ny, round(opts.rowRange(2))));
        rowsSelected = (min(rowA,rowB):max(rowA,rowB)).';
        firstPopulatedRow = rowsSelected(1);
        rowOrderSource = 'manual rowRange';
    else
        [rowsNearToFar, rowOrderSource] = camera_near_to_far_rows(Ny, opts, checkpoint);

        activeCountByRow = sum(mask, 2);
        minCells = opts.minPopulatedCellsPerRow;
        if isempty(minCells) || ~isfinite(minCells)
            % Avoid triggering on a single isolated populated cell while
            % still allowing narrow sensed regions.
            minCells = max(3, round(0.01 * Nx));
        end

        populated = activeCountByRow >= minCells;
        if ~any(populated)
            populated = activeCountByRow > 0;
        end

        firstIdxInOrder = find(populated(rowsNearToFar), 1, 'first');
        if isempty(firstIdxInOrder)
            warning('No populated WSE rows found in reference mask. Falling back to opts.rowRange.');
            rowA = max(1, min(Ny, round(opts.rowRange(1))));
            rowB = max(1, min(Ny, round(opts.rowRange(2))));
            rowsSelected = (min(rowA,rowB):max(rowA,rowB)).';
            firstPopulatedRow = rowsSelected(1);
            rowOrderSource = 'fallback manual rowRange';
        else
            firstPopulatedRow = rowsNearToFar(firstIdxInOrder);
            rowsAway = rowsNearToFar(firstIdxInOrder:end);

            dFromFirst_m = abs(yi(rowsAway) - yi(firstPopulatedRow));
            keep = dFromFirst_m <= opts.firstSensedLength_m;

            % Make sure at least one row is retained, and include the row
            % just beyond 10 m if the spacing skips over the boundary.
            if ~any(keep)
                keep(1) = true;
            end
            lastKeep = find(keep, 1, 'last');
            if lastKeep < numel(rowsAway) && dFromFirst_m(lastKeep) < opts.firstSensedLength_m
                keep(lastKeep + 1) = true;
            end

            rowsSelected = rowsAway(keep);
            rowsSelected = unique(rowsSelected(:), 'stable');
        end
    end

    if isempty(rowsSelected)
        error('No rows selected for stability assessment.');
    end

    nVals = yi(rowsSelected);
    targetN = make_axis(min(nVals), max(nVals), opts.profileSpacing_m);
    rowCoord = interp1(yi, 1:Ny, targetN, 'linear', 'extrap');

    activeCols = find(any(mask(rowsSelected,:), 1));
    if isempty(activeCols)
        activeCols = 1:Nx;
    end

    if isempty(opts.sampleSpacingAlongProfile_m)
        ds = median(abs(diff(xi)), 'omitnan');
    else
        ds = opts.sampleSpacingAlongProfile_m;
    end

    sLine = make_axis(min(xi(activeCols)), max(xi(activeCols)), ds);
    colCoord = interp1(xi, 1:Nx, sLine, 'linear', 'extrap');
    distance_m = sLine(:) - sLine(1);

    stabilityBand = struct;
    stabilityBand.mode = opts.stabilityAssessmentMode;
    stabilityBand.rows = rowsSelected(:);
    stabilityBand.rowOrderSource = rowOrderSource;
    stabilityBand.firstPopulatedRowIndex = firstPopulatedRow;
    stabilityBand.firstPopulatedRowCoord_m = yi(firstPopulatedRow);
    stabilityBand.startRowIndex = rowsSelected(1);
    stabilityBand.endRowIndex = rowsSelected(end);
    stabilityBand.startCoord_m = yi(rowsSelected(1));
    stabilityBand.endCoord_m = yi(rowsSelected(end));
    stabilityBand.length_m = abs(yi(rowsSelected(end)) - yi(rowsSelected(1)));
    stabilityBand.requestedLength_m = opts.firstSensedLength_m;
    stabilityBand.nRows = numel(rowsSelected);
    stabilityBand.nProfiles = numel(rowCoord);
    stabilityBand.minPopulatedCellsPerRow = opts.minPopulatedCellsPerRow;
end

function [rowsNearToFar, source] = camera_near_to_far_rows(Ny, opts, checkpoint)
    source = 'fallback descending rows';

    if strcmpi(char(opts.cameraRowOrder), 'ascending')
        rowsNearToFar = (1:Ny).';
        source = 'opts.cameraRowOrder=ascending';
        return
    elseif strcmpi(char(opts.cameraRowOrder), 'descending')
        rowsNearToFar = (Ny:-1:1).';
        source = 'opts.cameraRowOrder=descending';
        return
    end

    % In the solver checkpoint, y_order is the sweep order. The first entries
    % correspond to the near-camera side used during the downstream sweep.
    if isstruct(checkpoint) && isfield(checkpoint, 'y_order') && ~isempty(checkpoint.y_order)
        rowsNearToFar = round(checkpoint.y_order(:));
        rowsNearToFar = rowsNearToFar(isfinite(rowsNearToFar) & rowsNearToFar >= 1 & rowsNearToFar <= Ny);
        rowsNearToFar = unique(rowsNearToFar, 'stable');

        missingRows = setdiff((1:Ny).', rowsNearToFar(:), 'stable');
        rowsNearToFar = [rowsNearToFar(:); missingRows(:)];

        if ~isempty(rowsNearToFar)
            source = 'checkpoint.y_order';
            return
        end
    end

    rowsNearToFar = (Ny:-1:1).';
end

function x = make_axis(a, b, dx)
    if b < a
        tmp = a; a = b; b = tmp;
    end
    x = a:dx:b;
    if isempty(x) || x(end) < b
        x = [x, b];
    else
        x(end) = b;
    end
    x = x(:);
end

function z = sample_profile(W, rowCoord, colCoord)
    z = interp2(W, colCoord(:), rowCoord .* ones(size(colCoord(:))), 'linear', NaN);
end

function zDet = robust_linear_detrend(x, z)
    zDet = nan(size(z));
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
    zDet(good) = zDet(good) - mean(zDet(good), 'omitnan');
    zDet(~good) = NaN;
end

function a = robust_half_range(v, pLo, pHi)
    v = v(isfinite(v));
    if numel(v) < 5
        a = NaN;
    else
        a = 0.5 .* (prctile(v, pHi) - prctile(v, pLo));
    end
end

function [lambda, peakR] = autocorr_wavelength(x, z, minW, maxW, minR)
    lambda = NaN;
    peakR = NaN;
    good = isfinite(x) & isfinite(z);
    x = x(good);
    z = z(good);
    if numel(x) < 12
        return
    end
    [x, ord] = sort(x);
    z = z(ord);
    dx = median(diff(x), 'omitnan');
    if ~isfinite(dx) || dx <= 0
        return
    end
    z = robust_linear_detrend(x, z);
    good = isfinite(z);
    z = z(good);
    n = numel(z);
    if n < 12
        return
    end
    if ~isfinite(maxW)
        maxW = 0.75 * n * dx;
    end
    minLag = max(1, ceil(minW / dx));
    maxLag = min([floor(maxW / dx), floor(0.75*n), n-2]);
    if maxLag <= minLag
        return
    end
    ac = nan(maxLag,1);
    for lag = minLag:maxLag
        a = z(1:end-lag);
        b = z(1+lag:end);
        a = a - mean(a, 'omitnan');
        b = b - mean(b, 'omitnan');
        den = sqrt(sum(a.^2, 'omitnan') .* sum(b.^2, 'omitnan'));
        if den > 0
            ac(lag) = sum(a .* b, 'omitnan') ./ den;
        end
    end
    cand = [];
    for lag = minLag+1:maxLag-1
        if isfinite(ac(lag)) && ac(lag) >= minR && ac(lag) >= ac(lag-1) && ac(lag) >= ac(lag+1)
            cand(end+1) = lag; %#ok<AGROW>
        end
    end
    if isempty(cand)
        [peakR, lag] = max(ac(minLag:maxLag));
        lag = lag + minLag - 1;
        if ~isfinite(peakR) || peakR < minR
            return
        end
    else
        lag = cand(1);
        peakR = ac(lag);
    end
    delta = 0;
    if lag > 1 && lag < numel(ac) && all(isfinite(ac([lag-1 lag lag+1])))
        den = ac(lag-1) - 2*ac(lag) + ac(lag+1);
        if abs(den) > eps
            delta = max(-1, min(1, 0.5 * (ac(lag-1) - ac(lag+1)) ./ den));
        end
    end
    lambda = (lag + delta) * dx;
end


function R = pmusic_wavelength(x, z, minW, maxW, opts)
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

    good = isfinite(x) & isfinite(z);
    x = x(good);
    z = z(good);
    if numel(x) < 12
        R.status = 'too few samples';
        return
    end

    [x, ord] = sort(x(:));
    z = z(ord);
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

    zProcessed = robust_linear_detrend(xUniform, zUniform);
    goodZ = isfinite(zProcessed);
    if nnz(goodZ) < 12
        R.status = 'too few detrended samples';
        return
    end
    zProcessed = zProcessed(goodZ);
    zProcessed = zProcessed(:) - mean(zProcessed(:), 'omitnan');
    if sum(zProcessed.^2, 'omitnan') <= eps
        R.status = 'no variation';
        return
    end

    fsSpatial = 1 ./ dx; % samples per metre, so PMUSIC frequencies are cycles per metre
    modelOrder = max(1, round(opts.pmusicModelOrder));
    nfft = max(16, round(opts.pmusicNfft));

    if isempty(opts.pmusicWindowLength) || ~isfinite(opts.pmusicWindowLength)
        nwin = min(numel(zProcessed), max(modelOrder + 2, floor(numel(zProcessed)/2)));
    else
        nwin = min(numel(zProcessed), max(modelOrder + 2, round(opts.pmusicWindowLength)));
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

    if ~isfinite(maxW)
        maxW = max(minW, 0.75 .* (max(xUniform) - min(xUniform)));
    end

    fMin = 1 ./ maxW;
    fMax = 1 ./ minW;
    validPeak = isfinite(Pmusic) & isfinite(f) & f > 0 & f >= fMin & f <= fMax;

    if ~any(validPeak)
        R.status = 'no PMUSIC frequencies inside wavelength bounds';
        return
    end

    idxGood = find(validPeak);
    [peakPower, localIdx] = max(Pmusic(validPeak));
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

function v = winsorise(v, pLo, pHi)
    good = isfinite(v);
    if nnz(good) < 5
        return
    end
    q = prctile(v(good), [pLo pHi]);
    v(good) = max(q(1), min(q(2), v(good)));
end

function r = robust_rms(v)
    v = v(isfinite(v));
    if numel(v) < 5
        r = NaN;
        return
    end
    med = median(v, 'omitnan');
    madv = 1.4826 * median(abs(v - med), 'omitnan');
    if isfinite(madv) && madv > 0
        v = v(abs(v - med) <= 3.5*madv);
    end
    r = sqrt(mean(v.^2, 'omitnan'));
end

function [medv, q25, q75, n] = aggregate_profiles(V, goodMask, minProfiles)
    nMaps = size(V,1);
    medv = nan(nMaps,1);
    q25 = nan(nMaps,1);
    q75 = nan(nMaps,1);
    n = zeros(nMaps,1);
    for mm = 1:nMaps
        good = isfinite(V(mm,:)) & goodMask(mm,:);
        v = V(mm, good);
        if numel(v) >= minProfiles
            m = median(v, 'omitnan');
            madv = 1.4826 * median(abs(v-m), 'omitnan');
            if isfinite(madv) && madv > 0
                v = v(abs(v-m) <= 3.5*madv);
            end
        end
        n(mm) = numel(v);
        if ~isempty(v)
            medv(mm) = median(v, 'omitnan');
            q25(mm) = prctile(v, 25);
            q75(mm) = prctile(v, 75);
        end
    end
end

function M = metrics_to_map_series(checkpoint, nMaps)
    M = struct;
    M.fixedRepMedianAbsResidual = nan(nMaps,1);
    M.endRepCellMedianAbsResidual = nan(nMaps,1);
    M.localMedianRelativeImprovement = nan(nMaps,1);

    if ~isfield(checkpoint, 'metrics_static')
        return
    end
    ms = checkpoint.metrics_static;
    M.fixedRepMedianAbsResidual = hist_to_map(get_field_or_empty(ms, 'fixedRepPathMedianAbsResidualHist'), nMaps);
    M.endRepCellMedianAbsResidual = hist_to_map(get_field_or_empty(ms, 'endIterRepCellMedianAbsResidualHist'), nMaps);
    M.localMedianRelativeImprovement = hist_to_map(get_field_or_empty(ms, 'localMedianRelativeImprovementHist'), nMaps);
end

function x = get_field_or_empty(s, f)
    if isfield(s, f)
        x = s.(f);
    else
        x = [];
    end
end

function y = hist_to_map(h, nMaps)
    y = nan(nMaps,1);
    if isempty(h)
        return
    end
    h = h(:);
    n = min(numel(h), nMaps-1);
    y(2:n+1) = h(1:n);
end

function [selectedMap, method, info] = select_stop_map(Amed, Lmed, nA, nL, ampIqrFrac, wavIqrFrac, shapeChange, relUpdate, M, opts)
    nMaps = numel(Amed);
    info = struct;

    localMap = first_below(M.localMedianRelativeImprovement, opts.localImprovementCollapseThreshold, 4);
    fixedMap = first_near_min(M.fixedRepMedianAbsResidual, opts.residualNearMinimumFraction, 4);
    endMap = first_near_min(M.endRepCellMedianAbsResidual, opts.residualNearMinimumFraction, 4);

    info.localImprovementStopMap = localMap;
    info.fixedRep20PctMap = fixedMap;
    info.endRep20PctMap = endMap;
    info.requireForwardAmplitudePlateau = logical(opts.requireForwardAmplitudePlateau);
    info.noPlateauFallbackMode = char(opts.noPlateauFallbackMode);

    if isfinite(localMap)
        % Strong-signal branch.  Use residual-improvement collapse only as an
        % entrance gate.  The accepted map must also pass a forward amplitude
        % plateau test, so a declining relative map update alone cannot stop
        % the solver while the WSE amplitude is still growing.
        method = 'strong_signal_residual_collapse_plus_amplitude_plateau';
        selectedMap = NaN;
        fw = opts.strongForwardWindow;

        for s = localMap:nMaps
            sl = s:min(s+opts.stableWindow-1, nMaps);
            [qualityOk, qualityReason] = checkpoint_stop_quality_ok( ...
                Amed, Lmed, nA, nL, ampIqrFrac, wavIqrFrac, shapeChange, relUpdate, sl, 'strong', opts);
            if ~qualityOk
                continue
            end

            [plateauOk, forwardGrowth, forwardAbsGrowth_m, forwardTol_m] = amplitude_forward_plateau_ok( ...
                Amed, s, fw, opts.strongMaxForwardAmpGrowth, opts.strongMaxForwardAmpAbsGrowth_m);

            if plateauOk || ~opts.requireForwardAmplitudePlateau
                selectedMap = s;
                info.selectionReason = 'first strong-quality map satisfying amplitude plateau';
                info.qualityReason = qualityReason;
                info.forwardAmplitudeWindow = fw;
                info.forwardAmplitudeGrowth = forwardGrowth;
                info.forwardAmplitudeAbsGrowth_m = forwardAbsGrowth_m;
                info.forwardAmplitudeTolerance_m = forwardTol_m;
                return
            end
        end

        % If no plateau is visible within the available checkpoint sequence,
        % do not fall back to the first residual-collapse map.  That is the
        % failure mode that selected amplitudes around 0.02 m for synthetic
        % 0.06 m cases.  Instead, use the latest map that still satisfies the
        % amplitude/wavelength/shape quality checks.  This is conservative in
        % the sense that it only uses the final supported solution when there
        % is no evidence that amplitude has stabilised earlier.
        [fallbackMap, fallbackReason] = latest_quality_supported_map( ...
            Amed, Lmed, nA, nL, ampIqrFrac, wavIqrFrac, shapeChange, relUpdate, localMap, nMaps, 'strong', opts);
        if isfinite(fallbackMap)
            selectedMap = fallbackMap;
            method = 'strong_signal_no_plateau_latest_supported';
            info.selectionReason = fallbackReason;
            info.forwardAmplitudeWindow = fw;
            info.forwardAmplitudeGrowth = NaN;
            info.forwardAmplitudeAbsGrowth_m = NaN;
            info.forwardAmplitudeTolerance_m = NaN;
            return
        end
    end

    % Weak/low-signal branch.  The old version selected the first map where
    % relativeMapUpdate had decayed.  Because relativeMapUpdate = dA/A, it can
    % decay simply because A is growing, even when the absolute amplitude is
    % still increasing.  Therefore this branch now also requires a forward
    % amplitude plateau before accepting an early stop.
    method = 'weak_or_low_signal_update_decay_plus_amplitude_plateau';
    selectedMap = NaN;
    fw = opts.weakForwardWindow;
    for s = 2:nMaps
        sl = s:min(s+opts.stableWindow-1, nMaps);
        [qualityOk, qualityReason] = checkpoint_stop_quality_ok( ...
            Amed, Lmed, nA, nL, ampIqrFrac, wavIqrFrac, shapeChange, relUpdate, sl, 'weak', opts);
        if ~qualityOk
            continue
        end

        [plateauOk, forwardGrowth, forwardAbsGrowth_m, forwardTol_m] = amplitude_forward_plateau_ok( ...
            Amed, s, fw, opts.weakMaxForwardAmpGrowth, opts.weakMaxForwardAmpAbsGrowth_m);

        if plateauOk || ~opts.requireForwardAmplitudePlateau
            selectedMap = s;
            info.selectionReason = 'first weak-quality map satisfying amplitude plateau';
            info.qualityReason = qualityReason;
            info.forwardAmplitudeWindow = fw;
            info.forwardAmplitudeGrowth = forwardGrowth;
            info.forwardAmplitudeAbsGrowth_m = forwardAbsGrowth_m;
            info.forwardAmplitudeTolerance_m = forwardTol_m;
            return
        end
    end

    [fallbackMap, fallbackReason] = latest_quality_supported_map( ...
        Amed, Lmed, nA, nL, ampIqrFrac, wavIqrFrac, shapeChange, relUpdate, 2, nMaps, 'weak', opts);
    if isfinite(fallbackMap)
        selectedMap = fallbackMap;
        method = 'weak_or_low_signal_no_plateau_latest_supported';
        info.selectionReason = fallbackReason;
        info.forwardAmplitudeWindow = fw;
        info.forwardAmplitudeGrowth = NaN;
        info.forwardAmplitudeAbsGrowth_m = NaN;
        info.forwardAmplitudeTolerance_m = NaN;
        return
    end

    % Last-resort fallback: latest finite map with enough profile support.
    candidate = find(isfinite(Amed) & isfinite(Lmed) & nA >= opts.minProfiles & nL >= opts.minProfiles, 1, 'last');
    if isempty(candidate)
        error('No usable checkpoint map found.');
    end
    selectedMap = candidate;
    method = 'fallback_latest_supported_map';
    info.selectionReason = 'latest finite map with enough amplitude and wavelength profiles';
end

function [ok, reason] = checkpoint_stop_quality_ok(Amed, Lmed, nA, nL, ampIqrFrac, wavIqrFrac, shapeChange, relUpdate, sl, modeName, opts)
    ok = false;
    reason = 'not evaluated';

    sl = sl(:).';
    sl = sl(isfinite(sl));
    if numel(sl) < opts.stableWindow
        reason = 'insufficient forward stability window';
        return
    end

    if median(nA(sl), 'omitnan') < opts.minProfiles || median(nL(sl), 'omitnan') < opts.minProfiles
        reason = 'insufficient profile support';
        return
    end
    if any(~isfinite(Amed(sl))) || any(~isfinite(Lmed(sl)))
        reason = 'non-finite amplitude or wavelength in stability window';
        return
    end

    switch lower(char(modeName))
        case 'strong'
            ampIqrMax = opts.strongAmpIqrFracMax;
            wavIqrMax = opts.strongWavIqrFracMax;
            shapeTol = opts.strongShapeTol;
            checkRelUpdate = false;
        otherwise
            ampIqrMax = opts.weakAmpIqrFracMax;
            wavIqrMax = opts.weakWavIqrFracMax;
            shapeTol = opts.weakShapeTol;
            checkRelUpdate = true;
    end

    if median(ampIqrFrac(sl), 'omitnan') > ampIqrMax
        reason = 'amplitude IQR fraction too high';
        return
    end
    if median(wavIqrFrac(sl), 'omitnan') > wavIqrMax
        reason = 'wavelength IQR fraction too high';
        return
    end
    if median(shapeChange(sl), 'omitnan') > shapeTol
        reason = 'profile shape still changing';
        return
    end
    if checkRelUpdate && median(relUpdate(sl), 'omitnan') > opts.weakUpdateRelMax
        reason = 'relative map update has not decayed';
        return
    end

    Lwin = Lmed(sl);
    Lrange = prctile(Lwin, 90) - prctile(Lwin, 10);
    if Lrange > max(opts.wavStabAbs_m, opts.wavStabRel * median(Lwin, 'omitnan'))
        reason = 'wavelength not stable';
        return
    end

    ok = true;
    reason = 'quality checks passed';
end

function [ok, relativeGrowth, absoluteGrowth_m, tolerance_m] = amplitude_forward_plateau_ok(Amed, s, fw, relMax, absMax_m)
    ok = false;
    relativeGrowth = NaN;
    absoluteGrowth_m = NaN;
    tolerance_m = NaN;

    if ~isfinite(s) || ~isfinite(fw) || fw < 1
        return
    end
    s = round(s);
    fw = round(fw);
    if s < 1 || s+fw > numel(Amed)
        return
    end

    a0 = Amed(s);
    aFuture = Amed((s+1):(s+fw));
    if ~isfinite(a0) || any(~isfinite(aFuture))
        return
    end

    % Use the maximum future amplitude, not only the endpoint.  This catches
    % monotonic growth and also short overshoots that would be hidden by a
    % later decline.
    absoluteGrowth_m = max(aFuture) - a0;
    relativeGrowth = absoluteGrowth_m ./ max(abs(a0), eps);
    tolerance_m = max(absMax_m, relMax .* max(abs(a0), eps));
    ok = absoluteGrowth_m <= tolerance_m;
end

function [fallbackMap, reason] = latest_quality_supported_map(Amed, Lmed, nA, nL, ampIqrFrac, wavIqrFrac, shapeChange, relUpdate, startIdx, endIdx, modeName, opts)
    fallbackMap = NaN;
    reason = 'no quality-supported fallback map found';

    if ~isfinite(startIdx), startIdx = 1; end
    if ~isfinite(endIdx), endIdx = numel(Amed); end
    startIdx = max(1, round(startIdx));
    endIdx = min(numel(Amed), round(endIdx));

    if ~strcmpi(char(opts.noPlateauFallbackMode), 'latest_supported')
        reason = sprintf('no plateau found and noPlateauFallbackMode=%s', char(opts.noPlateauFallbackMode));
        return
    end

    for s = endIdx:-1:startIdx
        sl = max(startIdx, s-opts.stableWindow+1):s;
        [qualityOk, qualityReason] = checkpoint_stop_quality_ok( ...
            Amed, Lmed, nA, nL, ampIqrFrac, wavIqrFrac, shapeChange, relUpdate, sl, modeName, opts);
        if qualityOk
            fallbackMap = s;
            reason = sprintf('no forward amplitude plateau detected; using latest quality-supported map (%s)', qualityReason);
            return
        end
    end
end

function idx = first_below(y, thr, minIdx)
    idx = NaN;
    good = find(isfinite(y) & (y <= thr));
    good = good(good >= minIdx);
    if ~isempty(good)
        idx = good(1);
    end
end

function idx = first_near_min(y, frac, minIdx)
    idx = NaN;
    ys = movmedian(y, 5, 'omitnan');
    valid = find(isfinite(ys));
    valid = valid(valid >= minIdx);
    if isempty(valid)
        return
    end
    ymin = min(ys(valid));
    good = valid(ys(valid) <= ymin * (1 + frac));
    if ~isempty(good)
        idx = good(1);
    end
end
