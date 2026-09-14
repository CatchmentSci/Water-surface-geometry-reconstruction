function fig = KLT_plot_refined_wave_solution_three_panel_final( ...
        checkpointInput, refinedResult, mapTable, opts)
%KLT_PLOT_REFINED_WAVE_SOLUTION_THREE_PANEL_FINAL
% Create the final publication three-panel synthetic-wave figure.
%
% The figure contains:
%   (a) angular residuals;
%   (b) reconstructed WSE relative to the flat initial WSE; and
%   (c) imposed and reconstructed detrended WSE plus angular residuals
%       along the central transect of the first 10 m analysis band.
%
% The first 10 m analysis band begins at the first populated row on the
% camera side. Its middle available transect is plotted as the solid black
% line in (a) and (b) and supplies the data in (c). Panels (a) and (b) use
% the same spatial mask, limits and equal axis scaling. Masked cells have a
% white background. Panel (c) uses symmetric left and right limits so their
% zero levels align, and its angular residual is a solid orange line.
% Panel letters are positioned inside the axes to avoid PDF clipping.
%
% Usage with a checkpoint file:
%   fig = KLT_plot_refined_wave_solution_three_panel_final( ...
%       'syn_solver_inputs_0pt88_case_amp0.02_w1.52_ckpt.mat');
%
% Usage with selector outputs already in memory:
%   S = load('case_ckpt.mat','checkpoint');
%   [result,mapTable] = KLT_select_wave_checkpoint_refined(S.checkpoint);
%   fig = KLT_plot_refined_wave_solution_three_panel_final( ...
%       S.checkpoint,result,mapTable);
%
% Inputs:
%   checkpointInput  checkpoint struct or MAT-file containing 'checkpoint'
%   refinedResult    selector result; [] runs the selector internally
%   mapTable         selector map table; [] is accepted when the result
%                    contains selectedMapIndex
%   opts             optional settings structure
%
% Principal opts fields:
%   selectorOpts, transectProfileRange, profileToPlot, colormap,
%   figureWidthIn, figureHeightIn, saveFigure, outFigureFile,
%   usePhysicalAxesForMapPlots, mapXLimits, mapYLimits, clipMapXLimits,
%   clipMapYLimits, angularResidualCLim, relativeWseCLim, truthDEM,
%   syntheticTruthMode, syntheticTruthAmplitude_m and
%   syntheticTruthWavelength_m.
%
% Output and saving:
%   fig is the MATLAB figure handle. With saveFigure=true, an extensionless
%   outFigureFile writes both a 600-dpi PNG and a PDF. The default files are
%   example_syn.png and example_syn.pdf in the Synthetic analysis folder.

    if nargin < 1 || isempty(checkpointInput)
        if evalin('base', 'exist(''checkpoint'', ''var'')')
            checkpoint = evalin('base', 'checkpoint');
            checkpointFile = '';
        else
            error('Provide a checkpoint struct/file, or define checkpoint in the base workspace.');
        end
    elseif ischar(checkpointInput) || isstring(checkpointInput)
        checkpointFile = char(checkpointInput);
        S = load(checkpointFile, 'checkpoint');
        if ~isfield(S, 'checkpoint')
            error('The MAT file must contain a variable named checkpoint.');
        end
        checkpoint = S.checkpoint;
    elseif isstruct(checkpointInput)
        checkpoint = checkpointInput;
        checkpointFile = '';
    else
        error('checkpointInput must be a checkpoint struct, MAT filename, or empty.');
    end

    if nargin < 2, refinedResult = []; end
    if nargin < 3, mapTable = []; end
    if nargin < 4 || isempty(opts), opts = struct; end
    opts = set_plot_defaults(opts);

    if isempty(refinedResult)
        if exist('KLT_select_wave_checkpoint_refined', 'file') ~= 2
            error(['KLT_select_wave_checkpoint_refined.m was not found on the MATLAB path. ', ...
                   'Put it in the same folder or add its folder to the path.']);
        end
        [refinedResult, mapTable] = ...
            KLT_select_wave_checkpoint_refined(checkpoint, opts.selectorOpts);
    end

    if ~isfield(checkpoint, 'aa') || ~isfinite(checkpoint.aa)
        aa = 1;
    else
        aa = checkpoint.aa;
    end

    stopMap = resolve_stop_map(refinedResult, mapTable);

    if ~isfield(checkpoint, 'wse_map') || size(checkpoint.wse_map,1) < aa || size(checkpoint.wse_map,2) < stopMap
        error('checkpoint.wse_map{%d,%d} is not available.', aa, stopMap);
    end

    Waccepted = double(checkpoint.wse_map{aa, stopMap});
    [Ny, Nx] = size(Waccepted);

    if isfield(refinedResult, 'referenceMask') && isequal(size(refinedResult.referenceMask), size(Waccepted))
        referenceMask = logical(refinedResult.referenceMask);
    else
        referenceMask = isfinite(Waccepted);
    end

    angularResidualMap = get_residual_map( ...
        checkpoint, aa, stopMap, [Ny Nx]);
    angularDataMask = isfinite(angularResidualMap);
    displayMask = angularDataMask & referenceMask;

    % Relative-to-flat-WSE version for panel b.
    flatWseLevel = get_flat_wse_level(checkpoint, aa, referenceMask, Waccepted);
    WacceptedRelative = Waccepted - flatWseLevel;

    WacceptedDisplay = WacceptedRelative;
    if opts.maskAcceptedWseForDisplay
        WacceptedDisplay(~referenceMask) = NaN;
    end
    if opts.maskAcceptedWseToResidualCoverage
        WacceptedDisplay(~angularDataMask) = NaN;
    end

    if ~isfield(refinedResult, 'rowCoord') || ~isfield(refinedResult, 'colCoord') || ~isfield(refinedResult, 'distance_m')
        error('refinedResult must contain rowCoord, colCoord, and distance_m.');
    end
    rowCoord = refinedResult.rowCoord(:);
    colCoord = refinedResult.colCoord(:);
    distance_m = refinedResult.distance_m(:);
    nProfiles = numel(rowCoord);

    selectedProfileIndices = resolve_profile_indices(opts.transectProfileRange, nProfiles);
    if isempty(opts.profileToPlot)
        profileToPlot = selectedProfileIndices( ...
            ceil(numel(selectedProfileIndices) ./ 2));
    else
        profileToPlot = round(opts.profileToPlot);
    end
    if profileToPlot < 1 || profileToPlot > nProfiles
        error('profileToPlot=%d is outside available profile range 1:%d.', profileToPlot, nProfiles);
    end

    % Extract profile for panel c.
    nLine = min(numel(colCoord), numel(distance_m));
    colLine = colCoord(1:nLine);
    xRaw = distance_m(1:nLine);
    rowLine = rowCoord(profileToPlot) .* ones(nLine,1);

    estimatedWSE = interp2(Waccepted, colLine, rowLine, 'linear', NaN);
    angularResidualTransect = interp2(angularResidualMap, colLine, rowLine, 'linear', NaN);
    profileMask = interp2(double(displayMask), colLine, rowLine, 'nearest', 0) > 0.5;
    profileMask = profileMask(:) & isfinite(xRaw);

    if ~any(profileMask)
        warning('Central profile did not intersect finite residual/reference coverage. Using finite WSE coverage instead.');
        profileMask = isfinite(estimatedWSE) & isfinite(xRaw);
    end
    if ~any(profileMask)
        error('Profile %d has no finite WSE samples at the selected stop map.', profileToPlot);
    end

    idxFinite = find(profileMask);
    seg = idxFinite(1):idxFinite(end);
    distancePlot = xRaw(seg);
    distancePlot = distancePlot - distancePlot(1);

    estimatedWSE = estimatedWSE(seg);
    angularResidualTransect = angularResidualTransect(seg);
    mseg = profileMask(seg);
    estimatedWSE(~mseg) = NaN;
    angularResidualTransect(~mseg) = NaN;
    estimatedWSEDetrended = detrend_profile_for_amplitude(distancePlot, estimatedWSE, opts.detrendOrderForDisplay);

    truthWSEAmplitude = optional_truth_profile( ...
        checkpoint, checkpointFile, Waccepted, distance_m, colCoord, rowCoord(profileToPlot), seg, opts);

    [xMap, yMap, xLineAll, yLineAll, xLabelMap, yLabelMap] = make_map_axes(checkpoint, Nx, Ny, colCoord, rowCoord(profileToPlot), opts);
    xLinePlot = xLineAll(seg);
    yLinePlot = yLineAll(seg);

    selectedTransectY = rowCoord(selectedProfileIndices);
    if opts.usePhysicalAxesForMapPlots && isfield(checkpoint, 'yi') && numel(checkpoint.yi) == Ny
        selectedTransectY = interp1(1:Ny, yMap, selectedTransectY(:), 'linear', 'extrap');
    end

    [zoomXLim, zoomYLim] = compute_profile_region_limits(xLineAll, selectedTransectY, xLinePlot, yLinePlot, xMap, yMap, opts);

    fig = figure('Color','w', 'Units','inches', ...
        'Position',[0.5 0.5 opts.figureWidthIn opts.figureHeightIn], ...
        'PaperUnits','inches', 'PaperPosition',[0 0 opts.figureWidthIn opts.figureHeightIn], ...
        'PaperSize',[opts.figureWidthIn opts.figureHeightIn], 'InvertHardcopy','off');
    tl = tiledlayout(fig, 2, 2, 'TileSpacing','compact', 'Padding','compact');

    % Panel a
    ax1 = nexttile(tl, 1);
    h1 = imagesc(ax1, xMap, yMap, angularResidualMap);
    colormap(ax1, opts.colormap);
    set(h1, 'AlphaData', displayMask);
    set(ax1, 'YDir','normal', 'Color','w');
    hold(ax1, 'on');
    plot_transect_overlays(ax1, xLineAll, selectedTransectY, xLinePlot, yLinePlot, opts);
    cb1 = colorbar(ax1);
    cb1.Label.String = '$r_{\theta}$ ($^{\circ}$)';
    cb1.Label.Interpreter = 'latex';
    apply_zero_centred_clim(ax1, angularResidualMap(displayMask), opts.angularResidualCLim);
    apply_map_limits_and_equal(ax1, zoomXLim, zoomYLim, opts);
    xlabel(ax1, xLabelMap, 'Interpreter','latex');
    ylabel(ax1, yLabelMap, 'Interpreter','latex');
    title(ax1, '');
    style_map_axes(ax1, cb1, opts);
    add_panel_label(ax1, '(a)', opts.panelLabelPositionA, opts);

    % Panel b
    ax2 = nexttile(tl, 2);
    h2 = imagesc(ax2, xMap, yMap, WacceptedDisplay);
    colormap(ax2, opts.colormap);
    set(h2, 'AlphaData', isfinite(WacceptedDisplay));
    set(ax2, 'YDir','normal', 'Color','w');
    hold(ax2, 'on');
    plot_transect_overlays(ax2, xLineAll, selectedTransectY, xLinePlot, yLinePlot, opts);
    cb2 = colorbar(ax2);
    cb2.Label.String = '$\Delta z_{\mathrm{WSE}}$ (m)';
    cb2.Label.Interpreter = 'latex';
    apply_zero_centred_clim(ax2, WacceptedDisplay, opts.relativeWseCLim);
    apply_map_limits_and_equal(ax2, zoomXLim, zoomYLim, opts);
    xlabel(ax2, xLabelMap, 'Interpreter','latex');
    ylabel(ax2, yLabelMap, 'Interpreter','latex');
    title(ax2, '');
    style_map_axes(ax2, cb2, opts);
    add_panel_label(ax2, '(b)', opts.panelLabelPositionB, opts);


    % Force panels a) and b) to keep identical displayed extents.
    linkaxes([ax1, ax2], 'xy');
    % Use identical streamwise limits and ticks in panels (a) and (b).
    set([ax1 ax2], ...
        'XLim', [10 30], ...
        'XTick', 10:10:30, ...
        'YTick', -30:10:-10);
    ylim(ax2, ylim(ax1));

    % Panel c
    ax3 = nexttile(tl, [1 2]);
    hold(ax3, 'on');
    grid(ax3, 'on');
    set(ax3, 'FontSize',opts.tickFontSize, ...
        'LineWidth',opts.axesLineWidth, 'TickDir','out', ...
        'GridColor',opts.gridColor, 'GridAlpha',opts.gridAlpha, ...
        'Box','on', 'Layer','top');

    yyaxis(ax3, 'left');
    plotHandles = gobjects(0);
    plotLabels = {};

    if any(isfinite(truthWSEAmplitude))
        hTruth = plot(ax3, distancePlot, truthWSEAmplitude, '-', ...
            'LineWidth', 1.15, 'Color', [0 0.4470 0.7410], ...
            'DisplayName','Imposed geometry');
        plotHandles(end+1) = hTruth;
        plotLabels{end+1} = 'Imposed geometry';
    end

    hEstimated = plot(ax3, distancePlot, estimatedWSEDetrended, '-', ...
        'LineWidth', 1.25, 'Color', [0 0 0], ...
        'DisplayName','Reconstructed geometry');
    plotHandles(end+1) = hEstimated;
    plotLabels{end+1} = 'Reconstructed geometry';
    yline(ax3, 0, 'k:', 'LineWidth', 1.0, 'HandleVisibility','off');
    ylabel(ax3, '$z_{\mathrm{WSE}}^\prime$ (m)', ...
        'Interpreter','latex', 'FontSize',opts.axisLabelFontSize);

    leftVals = estimatedWSEDetrended(:);
    if any(isfinite(truthWSEAmplitude))
        leftVals = [leftVals; truthWSEAmplitude(:)];
    end
    leftYLim = symmetric_finite_ylim(leftVals);
    if ~isempty(leftYLim), ylim(ax3, leftYLim); end

    yyaxis(ax3, 'right');
    hAngular = plot(ax3, distancePlot, angularResidualTransect, '-', ...
        'LineWidth', 1.15, 'Color', [0.8500 0.3250 0.0980], ...
        'DisplayName','Angular residual');
    plotHandles(end+1) = hAngular;
    plotLabels{end+1} = 'Angular residual';
    ylabel(ax3, '$r_{\theta}$ ($^{\circ}$)', ...
        'Interpreter','latex', 'FontSize',opts.axisLabelFontSize);
    rightYLim = symmetric_finite_ylim(angularResidualTransect(:));
    if ~isempty(rightYLim), ylim(ax3, rightYLim); end
    ax3.YAxis(2).Color = [0.8500 0.3250 0.0980];

    yyaxis(ax3, 'left');
    ax3.YAxis(1).Color = [0 0 0];
    xlabel(ax3, '$s$ (m)', 'Interpreter','latex', ...
        'FontSize',opts.axisLabelFontSize);

    if isfinite(max(distancePlot)) && max(distancePlot) > 0
        xlim(ax3, [0 max(distancePlot)]);
    end
    lgd = legend(ax3, plotHandles, plotLabels, 'Location','southoutside', ...
        'Orientation','horizontal', 'Box','off', 'Interpreter','none');
    lgd.FontSize = opts.legendFontSize;
    title(ax3, '');
    add_panel_label(ax3, '(c)', opts.panelLabelPositionC, opts);

    if opts.saveFigure && ~isempty(opts.outFigureFile)
        [outDir,~,ext] = fileparts(opts.outFigureFile);
        if ~isempty(outDir) && ~isfolder(outDir)
            mkdir(outDir);
        end
        savedFigureFiles = {};
        if isempty(ext)
            pngFile = [opts.outFigureFile '.png'];
            pdfFile = [opts.outFigureFile '.pdf'];
            exportgraphics(fig, pngFile, 'Resolution', 600);
            exportgraphics(fig, pdfFile, ...
                'ContentType','vector', 'BackgroundColor','white');
            savedFigureFiles = {pngFile, pdfFile};
        elseif strcmpi(ext, '.pdf')
            exportgraphics(fig, opts.outFigureFile, ...
                'ContentType','vector', 'BackgroundColor','white');
            savedFigureFiles = {opts.outFigureFile};
        else
            exportgraphics(fig, opts.outFigureFile, 'Resolution', 600);
            savedFigureFiles = {opts.outFigureFile};
        end

        fprintf('\nFigure saved to:\n');
        for fileIndex = 1:numel(savedFigureFiles)
            [pathResolved, fileAttributes] = ...
                fileattrib(savedFigureFiles{fileIndex});
            if pathResolved
                fprintf('  %s\n', fileAttributes.Name);
            else
                fprintf('  %s\n', savedFigureFiles{fileIndex});
            end
        end
    end
end

%% ======================================================================
% Local helpers
%% ======================================================================
function opts = set_plot_defaults(opts)
    opts = default_field(opts, 'selectorOpts', struct());
    if isempty(opts.selectorOpts)
        opts.selectorOpts = struct;
    end
    % Only selector settings needed to reproduce the accepted map and
    % first-10-m transect geometry are specified here.
    opts.selectorOpts = default_field(opts.selectorOpts, 'stabilityAssessmentMode', 'first_populated_10m_from_camera');
    opts.selectorOpts = default_field(opts.selectorOpts, 'firstSensedLength_m', 10);
    opts.selectorOpts = default_field(opts.selectorOpts, 'cameraRowOrder', 'ascending');
    opts.selectorOpts = default_field(opts.selectorOpts, 'displaySummary', false);
    opts.selectorOpts = default_field(opts.selectorOpts, 'usePmusicWavelength', false);

    opts = default_field(opts, 'transectProfileRange', 'all');
    opts = default_field(opts, 'profileToPlot', []);
    opts = default_field(opts, 'colormap', local_coolwarm(256));
    opts = default_field(opts, 'figureWidthIn', 5.5);
    opts = default_field(opts, 'figureHeightIn', 4.4);
    opts = default_field(opts, 'saveFigure', true);
    opts = default_field(opts, 'outFigureFile', fullfile( ...
        fileparts(mfilename('fullpath')), 'output', 'Figure4'));
    opts = default_field(opts, 'maskAcceptedWseForDisplay', true);
    opts = default_field(opts, 'maskAcceptedWseToResidualCoverage', true);
    opts = default_field(opts, 'usePhysicalAxesForMapPlots', true);
    opts = default_field(opts, 'mapXLimits', []);
    opts = default_field(opts, 'mapYLimits', []);
    opts = default_field(opts, 'clipMapXLimits', false);
    opts = default_field(opts, 'clipMapYLimits', false);
    opts = default_field(opts, 'detrendOrderForDisplay', 1);
    opts = default_field(opts, 'transectColorAll', [0.4 0.4 0.4]);
    opts = default_field(opts, 'transectColorSelected', [0 0 0]);
    opts = default_field(opts, 'transectBandColor', [0.35 0.35 0.35]);
    opts = default_field(opts, 'transectBandAlpha', 0.08);
    opts = default_field(opts, 'transectBoundaryLineWidth', 0.55);
    opts = default_field(opts, 'transectZoomPaddingFraction', 0.08);
    opts = default_field(opts, 'axisLabelFontSize', 10);
    opts = default_field(opts, 'tickFontSize', 9);
    opts = default_field(opts, 'colorbarTickFontSize', 9);
    opts = default_field(opts, 'legendFontSize', 9);
    opts = default_field(opts, 'panelLabelFontSize', 10);
    opts = default_field(opts, 'axesLineWidth', 0.75);
    opts = default_field(opts, 'gridColor', [0.55 0.55 0.55]);
    opts = default_field(opts, 'gridAlpha', 0.14);
    opts = default_field(opts, 'angularResidualCLim', []);
    opts = default_field(opts, 'relativeWseCLim', []);
    opts = default_field(opts, 'syntheticTruthMode', 'filename_formula');
    opts = default_field(opts, 'syntheticTruthAmplitude_m', NaN);
    opts = default_field(opts, 'syntheticTruthWavelength_m', NaN);
    opts = default_field(opts, 'truthDEM', []);
    opts = default_field(opts, 'preferStoredTruth', true);
    opts = default_field(opts, 'panelLabelPositionA', [0.015 1.01423219727511]);
    opts = default_field(opts, 'panelLabelPositionB', [0.015 1.01423219727511]);
    opts = default_field(opts, 'panelLabelPositionC', [0.015 0.985]);
end

function flatWseLevel = get_flat_wse_level(checkpoint, aa, referenceMask, Waccepted)
    flatWseLevel = NaN;
    try
        if isfield(checkpoint, 'wse_map') && size(checkpoint.wse_map,1) >= aa && size(checkpoint.wse_map,2) >= 1
            W0 = checkpoint.wse_map{aa,1};
            if isnumeric(W0) && ~isempty(W0)
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

function profileIndices = resolve_profile_indices(profileRangeSetting, nProfiles)
    if ischar(profileRangeSetting) || isstring(profileRangeSetting)
        if strcmpi(char(profileRangeSetting), 'all')
            profileIndices = 1:nProfiles;
            return
        end
    end

    if isempty(profileRangeSetting)
        profileIndices = 1:nProfiles;
        return
    end

    profileRangeSetting = profileRangeSetting(:).';
    if numel(profileRangeSetting) == 2
        a = round(profileRangeSetting(1));
        b = round(profileRangeSetting(2));
        if b < a, tmp = a; a = b; b = tmp; end
        profileIndices = a:b;
    else
        profileIndices = unique(round(profileRangeSetting), 'stable');
    end
    profileIndices = profileIndices(isfinite(profileIndices) & profileIndices >= 1 & profileIndices <= nProfiles);
    if isempty(profileIndices)
        profileIndices = 1:nProfiles;
    end
end

function stopMap = resolve_stop_map(result, mapTable)
    if isfield(result, 'selectedMapRow') && isfinite(result.selectedMapRow)
        stopMap = round(result.selectedMapIndex);
        return
    end

    if ~isfield(result, 'selectedMapIndex') || ~isfinite(result.selectedMapIndex)
        error('result.selectedMapIndex is missing or non-finite.');
    end

    selectedValue = round(result.selectedMapIndex);
    stopMap = selectedValue;

    if ~isempty(mapTable) && istable(mapTable) && ismember('mapIndex', mapTable.Properties.VariableNames)
        exactRow = find(round(mapTable.mapIndex) == selectedValue, 1, 'first');
        if ~isempty(exactRow)
            stopMap = selectedValue;
            return
        end
        if selectedValue >= 1 && selectedValue <= height(mapTable)
            stopMap = round(mapTable.mapIndex(selectedValue));
            return
        end
    end
end

function M = get_residual_map(checkpoint, aa, stopMap, expectedSize)
    candidates = {};
    if isfield(checkpoint, 'metrics_static')
        ms = checkpoint.metrics_static;
        names = {'endIterRepResidualMapHist', 'fixedRepPathResidualMapHist', 'residualPostMapHist', 'residualMapHist'};
        for ii = 1:numel(names)
            if isfield(ms, names{ii})
                candidates{end+1} = ms.(names{ii}); %#ok<AGROW>
            end
        end
    end
    names2 = {'residualPostMapHist','residualMapHist','zi_post','zi','cell_median_val_init'};
    for ii = 1:numel(names2)
        if isfield(checkpoint, names2{ii})
            candidates{end+1} = checkpoint.(names2{ii}); %#ok<AGROW>
        end
    end

    % wse_map{...,1} is the initial surface. Residual history entry k is
    % evaluated on the surface after iteration k, i.e. wse_map{...,k+1}.
    expectedHistoryIndex = max(1, stopMap - 1);
    trialIdx = unique([expectedHistoryIndex stopMap stopMap+1 1], 'stable');
    for jj = 1:numel(trialIdx)
        for cc = 1:numel(candidates)
            M = extract_map_from_object( ...
                candidates{cc}, aa, trialIdx(jj), expectedSize);
            if ~isempty(M)
                if trialIdx(jj) ~= expectedHistoryIndex
                    warning(['Residual history index %d was unavailable; ', ...
                        'using fallback index %d.'], ...
                        expectedHistoryIndex, trialIdx(jj));
                end
                return
            end
        end
    end

    warning('No residual map found at the selected stop. Using finite WSE mask with NaN residual values.');
    M = NaN(expectedSize);
    M(isfinite(checkpoint.wse_map{aa, stopMap})) = 0;
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
            if ndims(obj) >= 2 && size(obj,1) >= aa && size(obj,2) >= mapIdx && isnumeric(obj{aa,mapIdx}) && isequal(size(obj{aa,mapIdx}), expectedSize)
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

function zDet = detrend_profile_for_amplitude(x, z, order)
    x = x(:); z = z(:);
    zDet = z;
    good = isfinite(x) & isfinite(z);
    if nnz(good) < 3
        return
    end
    if nargin < 3 || isempty(order) || ~isfinite(order)
        order = 1;
    end
    order = max(0, round(order));
    order = min(order, nnz(good)-1);
    if order == 0
        zDet(good) = z(good) - mean(z(good), 'omitnan');
        return
    end
    p = polyfit(x(good), z(good), order);
    zDet(good) = z(good) - polyval(p, x(good));
end

function [xMap, yMap, xLineAll, yLineAll, xLabelMap, yLabelMap] = make_map_axes(checkpoint, Nx, Ny, colCoord, rowValue, opts)
    xLabelMap = 'Column';
    yLabelMap = 'Row';
    if opts.usePhysicalAxesForMapPlots && isfield(checkpoint, 'xi') && numel(checkpoint.xi) == Nx && ...
            isfield(checkpoint, 'yi') && numel(checkpoint.yi) == Ny
        xMap = checkpoint.xi(:).';
        yMap = checkpoint.yi(:);
        xLineAll = interp1(1:Nx, xMap, colCoord(:), 'linear', 'extrap');
        yLineAll = interp1(1:Ny, yMap, rowValue(:).*ones(numel(colCoord),1), 'linear', 'extrap');
        xLabelMap = '$x$ (m)';
        yLabelMap = '$y$ (m)';
    else
        xMap = 1:Nx;
        yMap = 1:Ny;
        xLineAll = colCoord(:);
        yLineAll = rowValue(:).*ones(numel(colCoord),1);
    end
end

function [xLimOut, yLimOut] = compute_profile_region_limits(xLineAll, selectedTransectY, xLinePlot, yLinePlot, xMap, yMap, opts)
    xVals = [xLineAll(:); xLinePlot(:)];
    yVals = [selectedTransectY(:); yLinePlot(:)];

    xVals = xVals(isfinite(xVals));
    yVals = yVals(isfinite(yVals));
    if isempty(xVals)
        xVals = xMap(:);
    end
    if isempty(yVals)
        yVals = yMap(:);
    end

    xmin = min(xVals); xmax = max(xVals);
    ymin = min(yVals); ymax = max(yVals);

    xspan = xmax - xmin;
    yspan = ymax - ymin;

    if numel(xMap) > 1
        dx = median(diff(xMap), 'omitnan');
    else
        dx = 1;
    end
    if numel(yMap) > 1
        dy = median(diff(yMap), 'omitnan');
    else
        dy = 1;
    end

    if ~isfinite(xspan) || xspan <= 0, xspan = dx; end
    if ~isfinite(yspan) || yspan <= 0, yspan = dy; end

    xpad = max(opts.transectZoomPaddingFraction * xspan, 2*abs(dx));
    ypad = max(opts.transectZoomPaddingFraction * yspan, 2*abs(dy));

    xLimOut = [xmin - xpad, xmax + xpad];
    yLimOut = [ymin - ypad, ymax + ypad];

    xMapMin = min(xMap(:)); xMapMax = max(xMap(:));
    yMapMin = min(yMap(:)); yMapMax = max(yMap(:));
    xLimOut(1) = max(xLimOut(1), xMapMin);
    xLimOut(2) = min(xLimOut(2), xMapMax);
    yLimOut(1) = max(yLimOut(1), yMapMin);
    yLimOut(2) = min(yLimOut(2), yMapMax);
end

function plot_transect_overlays(ax, xLineAll, selectedTransectY, xLinePlot, yLinePlot, opts)
    finiteX = xLineAll(isfinite(xLineAll));
    finiteY = selectedTransectY(isfinite(selectedTransectY));

    if ~isempty(finiteX) && ~isempty(finiteY)
        xBounds = [min(finiteX), max(finiteX)];
        yBounds = [min(finiteY), max(finiteY)];
        patch(ax, ...
            [xBounds(1) xBounds(2) xBounds(2) xBounds(1)], ...
            [yBounds(1) yBounds(1) yBounds(2) yBounds(2)], ...
            opts.transectBandColor, ...
            'FaceAlpha',opts.transectBandAlpha, ...
            'EdgeColor','none', 'HandleVisibility','off');
        plot(ax, xBounds, [yBounds(1) yBounds(1)], '-', ...
            'Color',opts.transectColorAll, ...
            'LineWidth',opts.transectBoundaryLineWidth, ...
            'HandleVisibility','off');
        plot(ax, xBounds, [yBounds(2) yBounds(2)], '-', ...
            'Color',opts.transectColorAll, ...
            'LineWidth',opts.transectBoundaryLineWidth, ...
            'HandleVisibility','off');
    end
    plot(ax, xLinePlot, yLinePlot, '-', ...
        'Color',opts.transectColorSelected, 'LineWidth',1.2, ...
        'HandleVisibility','off');
end

function apply_zero_centred_clim(ax, M, requestedCLim)
    if nargin >= 3 && ~isempty(requestedCLim) && ...
            numel(requestedCLim) == 2 && all(isfinite(requestedCLim))
        requestedCLim = sort(double(requestedCLim(:).'));
        a = max(abs(requestedCLim));
        if a > 0
            clim(ax, [-a a]);
            return
        end
    end

    vals = M(isfinite(M));
    if isempty(vals), return; end
    a = max(abs(vals));
    if ~isfinite(a) || a == 0
        a = 1;
    end
    clim(ax, [-a a]);
end

function style_map_axes(ax, cb, opts)
    set(ax, 'FontSize',opts.tickFontSize, ...
        'LineWidth',opts.axesLineWidth, ...
        'TickDir','out', 'Box','on', 'Layer','top');
    ax.XLabel.FontSize = opts.axisLabelFontSize;
    ax.YLabel.FontSize = opts.axisLabelFontSize;
    cb.FontSize = opts.colorbarTickFontSize;
    cb.Label.FontSize = opts.axisLabelFontSize;
end

function apply_map_limits_and_equal(ax, zoomXLim, zoomYLim, opts)
    if opts.clipMapXLimits && ~isempty(opts.mapXLimits) && numel(opts.mapXLimits) == 2 && all(isfinite(opts.mapXLimits))
        xlim(ax, opts.mapXLimits);
    else
        xlim(ax, zoomXLim);
    end
    if opts.clipMapYLimits && ~isempty(opts.mapYLimits) && numel(opts.mapYLimits) == 2 && all(isfinite(opts.mapYLimits))
        ylim(ax, opts.mapYLimits);
    else
        ylim(ax, zoomYLim);
    end
    axis(ax, 'equal');
end

function yLim = symmetric_finite_ylim(vals)
    vals = vals(:);
    vals = vals(isfinite(vals));
    if isempty(vals)
        yLim = [];
        return
    end
    a = max(abs(vals));
    if ~isfinite(a) || a == 0
        a = 1;
    end
    yLim = [-a a];
end

function add_panel_label(ax, txt, position, opts)
    text(ax, position(1), position(2), txt, ...
        'Units','normalized', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','top', ...
        'FontWeight','bold', ...
        'FontSize',opts.panelLabelFontSize, ...
        'Color','k', ...
        'Interpreter','none', ...
        'Clipping','on');
end

function truthAmp = optional_truth_profile(checkpoint, checkpointFile, Waccepted, distance_m, colCoord, rowProfile, profileSegmentIdx, opts)
    truthAmp = NaN(numel(profileSegmentIdx),1);
    mode = lower(char(opts.syntheticTruthMode));
    truthDEM = [];

    if opts.preferStoredTruth
        truthDEM = find_stored_truth_dem( ...
            checkpoint, opts, size(Waccepted));
    end

    if isempty(truthDEM) && (strcmp(mode, 'manual_formula') || ...
            strcmp(mode, 'filename_formula') || strcmp(mode, 'formula'))
        amp = opts.syntheticTruthAmplitude_m;
        wln = opts.syntheticTruthWavelength_m;
        if strcmp(mode, 'filename_formula') && (~isfinite(amp) || ~isfinite(wln))
            [amp2, wln2] = parse_amp_wavelength_from_name(checkpointFile);
            if ~isfinite(amp), amp = amp2; end
            if ~isfinite(wln), wln = wln2; end
        end
        if isfinite(amp) && isfinite(wln) && wln > 0 && isfield(checkpoint, 'xi') && isfield(checkpoint, 'yi')
            [Sg, ~] = meshgrid(checkpoint.xi(:).', checkpoint.yi(:));
            truthDEM = amp .* sin(2*pi*Sg ./ wln);
        end
    end

    if isempty(truthDEM)
        return
    end

    nLine = min(numel(colCoord), numel(distance_m));
    rowLine = rowProfile .* ones(nLine,1);
    truthWSE = interp2(truthDEM, colCoord(1:nLine), rowLine, 'linear', NaN);
    truthWSE = truthWSE(profileSegmentIdx);
    dist = distance_m(profileSegmentIdx);
    dist = dist - dist(1);
    truthAmp = detrend_profile_for_amplitude(dist, truthWSE, 1);
end

function truthDEM = find_stored_truth_dem(checkpoint, opts, expectedSize)
    truthDEM = [];

    if ~isempty(opts.truthDEM) && isnumeric(opts.truthDEM) && ...
            isequal(size(opts.truthDEM), expectedSize)
        truthDEM = double(opts.truthDEM);
        return
    end

    candidateFields = { ...
        'truthDEM', ...
        'truthWSE', ...
        'syntheticTruthDEM', ...
        'syntheticTruthWSE', ...
        'wseTruth', ...
        'WSE_truth'};

    for ii = 1:numel(candidateFields)
        name = candidateFields{ii};
        if isfield(checkpoint, name)
            value = checkpoint.(name);
            if isnumeric(value) && isequal(size(value), expectedSize)
                truthDEM = double(value);
                return
            end
        end
    end
end

function [amp, wln] = parse_amp_wavelength_from_name(fileName)
    amp = NaN;
    wln = NaN;
    if isempty(fileName)
        return
    end
    [~, name, ext] = fileparts(fileName);
    s = [name ext];
    tok = regexp(s, 'amp([0-9]+(?:\.[0-9]+)?)_w([0-9]+(?:\.[0-9]+)?)', 'tokens', 'once');
    if ~isempty(tok)
        amp = str2double(tok{1});
        wln = str2double(tok{2});
    end
end

function s = default_field(s, f, v)
    if ~isfield(s, f) || isempty(s.(f))
        s.(f) = v;
    end
end

function map = local_coolwarm(m)
    if nargin < 1 || isempty(m)
        m = 256;
    end
    if m <= 0
        map = zeros(0,3);
        return
    end
    % Compact coolwarm anchor set, interpolated to requested size.
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
    x = linspace(0,1,size(anchors,1));
    xi = linspace(0,1,m);
    map = interp1(x, anchors, xi, 'linear');
    map = max(0, min(1, map));
end
