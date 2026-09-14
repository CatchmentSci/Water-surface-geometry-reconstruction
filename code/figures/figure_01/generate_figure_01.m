function outputs = generate_figure_01(dataRoot, outputDir)
%GENERATE_FIGURE_01 Reproduce Figure 1 from the archived R9 solver data.
% -------------------------------------------------------------------------
% Plot the WSE reconstruction focus region and flow vectors in image space.
%
% This version:
%   1) Builds the WSE mask from cells actually adjusted by the solver using
%      checkpoint.dzMapHist.
%   2) Projects the unsimplified adjusted-cell boundary into image space
%      using the default/background Z value, e.g. 4999.19.
%   3) Rasterises the projected polygon into an image-space mask using
%      poly2mask.
%   4) Selects vectors by testing their image-space A/B/midpoint pixel
%      positions against that raster mask. This avoids slow inverse-projection
%      of millions of vectors.
%   5) Colours plotted vectors by pixel-space angular deviation from the
%      idealised flow direction defined by app.pts.
%   6) Sets figure width to 17.8 cm and all fonts to Arial, size 12.
% -------------------------------------------------------------------------

%% ------------------------------------------------------------------------
% Inputs and repository paths
% -------------------------------------------------------------------------

arguments
    dataRoot (1,1) string
    outputDir (1,1) string = string(fullfile(fileparts(mfilename('fullpath')), 'output'))
end

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(fileparts(scriptDir)));
dependencyDir = fullfile(repoRoot, 'code', 'dependencies');
addpath(dependencyDir);

inputFile = fullfile(dataRoot, 'videos', 'inputs', ...
    'devon_dart20180315_08180_solver_inputs.mat');
checkpointFile = fullfile(dataRoot, 'videos', 'outputs', ...
    'devon_dart20180315_08180_checkpoint.mat');

if ~isfile(inputFile)
    error('Figure 1 solver input file not found: %s', inputFile);
end
if ~isfile(checkpointFile)
    error('Figure 1 checkpoint file not found: %s', checkpointFile);
end
if ~isfolder(outputDir)
    mkdir(outputDir);
end

inputData = load(inputFile, 'app_in', 'camA_fullmodel', 'xyzA_wse', 'xyzB_wse');
checkpointData = load(checkpointFile, 'checkpoint');

requiredInputs = {'app_in', 'camA_fullmodel', 'xyzA_wse', 'xyzB_wse'};
for k = 1:numel(requiredInputs)
    if ~isfield(inputData, requiredInputs{k})
        error('Variable %s is missing from %s.', requiredInputs{k}, inputFile);
    end
end
if ~isfield(checkpointData, 'checkpoint')
    error('Variable checkpoint is missing from %s.', checkpointFile);
end

app = inputData.app_in;
camA = camera(inputData.camA_fullmodel);
app.camA = camA;
xyzA_wse = inputData.xyzA_wse;
xyzB_wse = inputData.xyzB_wse;
checkpoint = checkpointData.checkpoint;

aa = 1;
vectorCellIdx = 1;

% Use [] to use final available WSE map
mapIdxUser = [];

% Default/background Z used for polygon projection
defaultProjectionZ = 4999.19;

% Adjusted-cell mask settings
updateTol = 0;              % non-zero dz threshold
minBlobSizeCells = 20;
fillMaskHoles = false;      % false = only actually adjusted cells
keepLargestBlobOnly = true; % true = one main polygon

% Vector inclusion rule:
%   'midpoint' = vector midpoint inside projected adjusted WSE mask
%   'A'        = vector start point inside projected adjusted WSE mask
%   'B'        = vector end point inside projected adjusted WSE mask
%   'either'   = start or end point inside projected adjusted WSE mask
%   'both'     = start and end point inside projected adjusted WSE mask
vectorInsideMode = 'midpoint';

% Plotting control
maxVectorsToPlot = 200000;   % use Inf to plot all selected vectors

% Pixel-space angle colouring
colourVectorsByPixelAngle = true;
wrapAngleTo180 = true;

% Set to 18.5 if you want to subtract 18.5 deg from plotted angle deviations.
% Set to 0 for pure deviation from app.pts.
pixelAngleCorrectionDeg = 0;

% Colour limits in degrees. Use [] for automatic symmetric limits.
angleCLim = [-40 40];

% Repository-owned approximation of the cool-warm diverging colour map.
nAngleColours = 256;
angleColormap = figure_01_colormap(nAngleColours);

polygonColor = [1 0 0];
vectorLineWidth = 1;
polygonLineWidth = 1;

% Figure / font formatting
plotFontName = 'Arial';
plotFontSize = 10;

figureWidthCm = 18;
figureHeightCm = [];   % [] = calculate from image aspect ratio
figureLeftCm = 2;
figureBottomCm = 2;

% Output
outputFile = fullfile(outputDir, ...
    'projected_wse_adjusted_region_with_pixel_angle_vectors_fast.mat');
figureFile = fullfile(outputDir, 'Figure1.png');

%% ------------------------------------------------------------------------
% Select WSE map
% -------------------------------------------------------------------------

if ~isfield(checkpoint, 'wse_map')
    error('checkpoint.wse_map is missing.');
end

if isempty(mapIdxUser)
    nonEmptyMaps = find(~cellfun(@isempty, checkpoint.wse_map(aa,:)));

    if isempty(nonEmptyMaps)
        error('No non-empty WSE maps found in checkpoint.wse_map for aa = %d.', aa);
    end

    mapIdx = nonEmptyMaps(end);
else
    mapIdx = mapIdxUser;
end

W = checkpoint.wse_map{aa, mapIdx};

if isempty(W)
    error('checkpoint.wse_map{%d,%d} is empty.', aa, mapIdx);
end

fprintf('Using checkpoint.wse_map{%d,%d}\n', aa, mapIdx);
fprintf('WSE map size: %d rows x %d columns\n', size(W,1), size(W,2));

%% ------------------------------------------------------------------------
% Check required checkpoint fields
% -------------------------------------------------------------------------

requiredFields = {'xi', 'yi', 'origin', 'e_s', 'e_n'};

for k = 1:numel(requiredFields)
    f = requiredFields{k};

    if ~isfield(checkpoint, f)
        error('checkpoint.%s is missing.', f);
    end
end

Ny = size(W, 1);
Nx = size(W, 2);

if numel(checkpoint.xi) ~= Nx
    error('numel(checkpoint.xi) does not match number of WSE columns.');
end

if numel(checkpoint.yi) ~= Ny
    error('numel(checkpoint.yi) does not match number of WSE rows.');
end

%% ------------------------------------------------------------------------
% Build adjusted-cell WSE mask from dzMapHist
% -------------------------------------------------------------------------

dataMask = false(size(W));

if isfield(checkpoint, 'dzMapHist') && ~isempty(checkpoint.dzMapHist)
    dzHist = checkpoint.dzMapHist;

elseif isfield(checkpoint, 'metrics_static') && ...
        isfield(checkpoint.metrics_static, 'dzMapHist') && ...
        ~isempty(checkpoint.metrics_static.dzMapHist)
    dzHist = checkpoint.metrics_static.dzMapHist;

elseif isfield(checkpoint, 'metrics') && ...
        isfield(checkpoint.metrics, 'dzMapHist') && ...
        ~isempty(checkpoint.metrics.dzMapHist)
    dzHist = checkpoint.metrics.dzMapHist;

else
    dzHist = {};
end

if ~isempty(dzHist)

    maxDzIter = min(numel(dzHist), max(mapIdx - 1, 0));

    if maxDzIter < 1
        error('Selected mapIdx=%d has no corresponding dzMapHist updates.', mapIdx);
    end

    for k = 1:maxDzIter
        dz = dzHist{k};

        if isempty(dz)
            continue
        end

        if ~isequal(size(dz), size(W))
            warning('Skipping dzMapHist{%d}: size does not match WSE map.', k);
            continue
        end

        dataMask = dataMask | (isfinite(dz) & abs(dz) > updateTol);
    end

    fprintf('Adjusted-cell mask built from dzMapHist iterations 1:%d\n', maxDzIter);

else

    warning(['No dzMapHist found. Falling back to WSE difference between ', ...
             'initial and selected WSE map.']);

    W0 = checkpoint.wse_map{aa,1};

    if isempty(W0) || ~isequal(size(W0), size(W))
        error('Cannot build fallback mask: checkpoint.wse_map{%d,1} is missing or wrong size.', aa);
    end

    dataMask = isfinite(W) & isfinite(W0) & abs(W - W0) > updateTol;
end

fprintf('Raw adjusted-cell mask cells = %d / %d (%.3f%%)\n', ...
    nnz(dataMask), numel(dataMask), 100 * nnz(dataMask) / numel(dataMask));

if minBlobSizeCells > 0
    dataMask = bwareaopen(dataMask, minBlobSizeCells);
end

if fillMaskHoles
    dataMask = imfill(dataMask, 'holes');
end

if keepLargestBlobOnly && any(dataMask(:))
    dataMask = bwareafilt(dataMask, 1);
end

if ~any(dataMask(:))
    error('Adjusted-cell dataMask is empty after cleaning.');
end

[rowIdx, colIdx] = find(dataMask);

fprintf('Cleaned adjusted-cell mask cells = %d / %d (%.3f%%)\n', ...
    nnz(dataMask), numel(dataMask), 100 * nnz(dataMask) / numel(dataMask));
fprintf('Adjusted-cell row range: %d:%d\n', min(rowIdx), max(rowIdx));
fprintf('Adjusted-cell col range: %d:%d\n', min(colIdx), max(colIdx));

%% ------------------------------------------------------------------------
% Extract actual unsimplified boundary/boundaries from adjusted-cell mask
% -------------------------------------------------------------------------

B = bwboundaries(dataMask, 8, 'noholes');

if isempty(B)
    error('No boundary found for adjusted-cell dataMask.');
end

fprintf('Number of adjusted-cell boundaries: %d\n', numel(B));

%% ------------------------------------------------------------------------
% Build image-space polygon boundary/boundaries using default Z
% -------------------------------------------------------------------------

zPlane = defaultProjectionZ;

origin = checkpoint.origin;
e_s = checkpoint.e_s;
e_n = checkpoint.e_n;

imageBoundaryUVCell = cell(numel(B), 1);
imageBoundaryFiniteMaskCell = cell(numel(B), 1);
worldBoundaryXYZCell = cell(numel(B), 1);
rotatedBoundarySNCell = cell(numel(B), 1);
wseBoundaryColRowCell = cell(numel(B), 1);

polyUVPlotAll = [];
worldBoundaryXYZAll = [];
rotatedBoundarySNAll = [];
wseBoundaryColRowAll = [];

fprintf('Using default/background projection height Z = %.6f m\n', zPlane);

for ib = 1:numel(B)

    boundaryRC = B{ib};   % [row, col]

    rowB = boundaryRC(:,1);
    colB = boundaryRC(:,2);

    if colB(1) ~= colB(end) || rowB(1) ~= rowB(end)
        colB(end+1) = colB(1);
        rowB(end+1) = rowB(1);
    end

    sB = interp1(1:numel(checkpoint.xi), checkpoint.xi, colB, 'linear', 'extrap');
    nB = interp1(1:numel(checkpoint.yi), checkpoint.yi, rowB, 'linear', 'extrap');

    XB = origin(1) + sB(:).*e_s(1) + nB(:).*e_n(1);
    YB = origin(2) + sB(:).*e_s(2) + nB(:).*e_n(2);
    ZB = zPlane * ones(numel(sB), 1);

    polyXYZ = [XB(:), YB(:), ZB(:)];

    try
        [polyUV, ~, inFrame] = camA.project(polyXYZ); %#ok<ASGLU>
    catch
        [polyUV, ~, inFrame] = project(camA, polyXYZ); %#ok<ASGLU>
    end

    finiteUV = all(isfinite(polyUV), 2);

    imageBoundaryUVCell{ib} = polyUV;
    imageBoundaryFiniteMaskCell{ib} = finiteUV;
    worldBoundaryXYZCell{ib} = polyXYZ;
    rotatedBoundarySNCell{ib} = [sB(:), nB(:)];
    wseBoundaryColRowCell{ib} = [colB(:), rowB(:)];

    polyUVPlot = polyUV;
    polyUVPlot(~finiteUV,:) = NaN;

    polyUVPlotAll = [polyUVPlotAll; polyUVPlot; NaN NaN]; %#ok<AGROW>
    worldBoundaryXYZAll = [worldBoundaryXYZAll; polyXYZ; NaN NaN NaN]; %#ok<AGROW>
    rotatedBoundarySNAll = [rotatedBoundarySNAll; [sB(:), nB(:)]; NaN NaN]; %#ok<AGROW>
    wseBoundaryColRowAll = [wseBoundaryColRowAll; [colB(:), rowB(:)]; NaN NaN]; %#ok<AGROW>

    fprintf('Boundary %d: vertices=%d, finite projected=%d/%d\n', ...
        ib, numel(colB), nnz(finiteUV), numel(finiteUV));
end

%% ------------------------------------------------------------------------
% Get image-space flow-vector endpoints
% -------------------------------------------------------------------------

uvA = local_get_wse_vector_cell(xyzA_wse, aa, vectorCellIdx, 'xyzA_wse');
uvB = local_get_wse_vector_cell(xyzB_wse, aa, vectorCellIdx, 'xyzB_wse');

uvA = uvA(:,1:2);
uvB = uvB(:,1:2);

if size(uvA,1) ~= size(uvB,1)
    error('uvA and uvB must have the same number of rows.');
end

finiteVec = all(isfinite(uvA), 2) & all(isfinite(uvB), 2);

xA = uvA(:,1);
yA = uvA(:,2);
xB = uvB(:,1);
yB = uvB(:,2);

xMid = 0.5 * (xA + xB);
yMid = 0.5 * (yA + yB);

fprintf('Total vectors: %d\n', size(uvA,1));
fprintf('Finite vectors: %d\n', nnz(finiteVec));

%% ------------------------------------------------------------------------
% Select vectors using image-space raster mask
% -------------------------------------------------------------------------
% This is the fast replacement for inverse-projecting vector midpoints.
% The adjusted WSE region has already been projected into image space.

fprintf('Selecting vectors using image-space raster mask from projected WSE polygon.\n');

if isfield(app, 'firstFrame') && ~isempty(app.firstFrame)
    imageRows = size(app.firstFrame, 1);
    imageCols = size(app.firstFrame, 2);
else
    allU = [uvA(:,1); uvB(:,1); polyUVPlotAll(:,1)];
    allV = [uvA(:,2); uvB(:,2); polyUVPlotAll(:,2)];

    imageCols = ceil(max(allU(isfinite(allU))));
    imageRows = ceil(max(allV(isfinite(allV))));

    warning('app.firstFrame missing; inferred image size as %d rows x %d cols.', ...
        imageRows, imageCols);
end

imageRegionMask = false(imageRows, imageCols);

for ib = 1:numel(imageBoundaryUVCell)

    uvPoly = imageBoundaryUVCell{ib};

    if isempty(uvPoly)
        continue
    end

    goodPoly = all(isfinite(uvPoly), 2);

    if nnz(goodPoly) < 3
        continue
    end

    xPoly = uvPoly(goodPoly, 1);
    yPoly = uvPoly(goodPoly, 2);

    imageRegionMask = imageRegionMask | poly2mask(xPoly, yPoly, imageRows, imageCols);
end

fprintf('Image-space polygon mask pixels = %d / %d (%.3f%%)\n', ...
    nnz(imageRegionMask), numel(imageRegionMask), ...
    100 * nnz(imageRegionMask) / numel(imageRegionMask));

switch lower(vectorInsideMode)

    case 'midpoint'
        inRegion = local_image_points_inside_mask( ...
            xMid, yMid, finiteVec, imageRegionMask);

    case 'a'
        inRegion = local_image_points_inside_mask( ...
            xA, yA, finiteVec, imageRegionMask);

    case 'b'
        inRegion = local_image_points_inside_mask( ...
            xB, yB, finiteVec, imageRegionMask);

    case 'either'
        inA = local_image_points_inside_mask( ...
            xA, yA, finiteVec, imageRegionMask);

        inB = local_image_points_inside_mask( ...
            xB, yB, finiteVec, imageRegionMask);

        inRegion = finiteVec & (inA | inB);

    case 'both'
        inA = local_image_points_inside_mask( ...
            xA, yA, finiteVec, imageRegionMask);

        inB = local_image_points_inside_mask( ...
            xB, yB, finiteVec, imageRegionMask);

        inRegion = finiteVec & inA & inB;

    otherwise
        error('Unknown vectorInsideMode: %s', vectorInsideMode);
end

idxVectorsInRegion = find(inRegion);

fprintf('Vectors inside projected adjusted-cell polygon using mode "%s": %d\n', ...
    vectorInsideMode, numel(idxVectorsInRegion));

fprintf('Selected vector fraction = %.3f%%\n', ...
    100 * numel(idxVectorsInRegion) / max(nnz(finiteVec), 1));

%% ------------------------------------------------------------------------
% Subsample vectors for plotting
% -------------------------------------------------------------------------

idxVectorsToPlot = idxVectorsInRegion;

if isfinite(maxVectorsToPlot) && numel(idxVectorsToPlot) > maxVectorsToPlot
    keepSub = round(linspace(1, numel(idxVectorsToPlot), maxVectorsToPlot));
    idxVectorsToPlot = idxVectorsToPlot(keepSub);

    fprintf('Subsampled vectors for plotting: %d of %d\n', ...
        numel(idxVectorsToPlot), numel(idxVectorsInRegion));
else
    fprintf('Plotting all selected vectors: %d\n', numel(idxVectorsToPlot));
end

%% ------------------------------------------------------------------------
% Compute pixel-space angular deviation for plotted vectors
% -------------------------------------------------------------------------

pixelAngleDeviationPlot = [];

if ~isempty(idxVectorsToPlot)

    if ~isfield(app, 'pts') || isempty(app.pts) || ...
            size(app.pts,1) < 2 || size(app.pts,2) < 2
        error('app.pts must contain at least two [x;y] points defining the ideal pixel-space flow direction.');
    end

    startPix = [app.pts(1,1), app.pts(2,1)];
    endPix   = [app.pts(1,2), app.pts(2,2)];

    idealVecPix = endPix - startPix;

    if ~all(isfinite(idealVecPix)) || hypot(idealVecPix(1), idealVecPix(2)) == 0
        error('Ideal pixel-space flow direction from app.pts is invalid.');
    end

    idealAnglePix = atan2d(idealVecPix(2), idealVecPix(1));

    dxPix = xB(idxVectorsToPlot) - xA(idxVectorsToPlot);
    dyPix = yB(idxVectorsToPlot) - yA(idxVectorsToPlot);

    obsAnglePix = atan2d(dyPix, dxPix);

    pixelAngleDeviationPlot = idealAnglePix - obsAnglePix;

    if wrapAngleTo180
        pixelAngleDeviationPlot = local_wrap180(pixelAngleDeviationPlot);
    end

    pixelAngleDeviationPlot = pixelAngleDeviationPlot - pixelAngleCorrectionDeg;

    fprintf('Pixel-space ideal flow angle from app.pts = %+0.3f deg\n', idealAnglePix);
    fprintf('Computed pixel-space angular deviations for %d plotted vectors.\n', ...
        numel(idxVectorsToPlot));
end

%% ------------------------------------------------------------------------
% Plot image with adjusted-cell polygon and selected vectors
% -------------------------------------------------------------------------

% Work out figure height from image aspect ratio
if isempty(figureHeightCm)
    if exist('imageRows', 'var') && exist('imageCols', 'var') && ...
            isfinite(imageRows) && isfinite(imageCols) && imageCols > 0
        figureHeightCm = figureWidthCm * imageRows / imageCols;
    else
        figureHeightCm = figureWidthCm * 0.65;
    end
end

figH = figure( ...
    'Color', 'w', ...
    'Units', 'centimeters', ...
    'Position', [figureLeftCm, figureBottomCm, figureWidthCm, figureHeightCm], ...
    'PaperUnits', 'centimeters', ...
    'PaperPosition', [0, 0, figureWidthCm, figureHeightCm], ...
    'PaperSize', [figureWidthCm, figureHeightCm]);

if isfield(app, 'firstFrame') && ~isempty(app.firstFrame)

    bgFrame = app.firstFrame;

    % Convert grayscale / indexed image to RGB so residual colormap does not
    % recolour the background.
    if ndims(bgFrame) == 2
        bgFrame = mat2gray(bgFrame);
        bgFrame = repmat(bgFrame, 1, 1, 3);
    elseif ndims(bgFrame) == 3 && size(bgFrame,3) == 1
        bgFrame = mat2gray(bgFrame(:,:,1));
        bgFrame = repmat(bgFrame, 1, 1, 3);
    elseif ndims(bgFrame) == 3 && size(bgFrame,3) == 3
        % Already RGB. Leave as-is.
    else
        error('app.firstFrame has unsupported dimensions.');
    end

    imshow(bgFrame);
    hold on;

else

    hold on;
    axis equal;
    set(gca, 'YDir', 'reverse');
    grid on;
    xlabel('Image column / u');
    ylabel('Image row / v');

end

if ~isempty(idxVectorsToPlot)

    if colourVectorsByPixelAngle && ~isempty(pixelAngleDeviationPlot)

        goodAngle = isfinite(pixelAngleDeviationPlot);

        idxGoodPlot = idxVectorsToPlot(goodAngle);
        angleGood = pixelAngleDeviationPlot(goodAngle);

        if ~isempty(idxGoodPlot)

            cmap = angleColormap;
            nColours = size(cmap, 1);

            if isempty(angleCLim)
                cmax = prctile(abs(angleGood), 98);

                if ~isfinite(cmax) || cmax == 0
                    cmax = 1;
                end

                climUse = [-cmax cmax];
            else
                climUse = angleCLim;
            end

            anglePlot = angleGood;
            anglePlot(anglePlot < climUse(1)) = climUse(1);
            anglePlot(anglePlot > climUse(2)) = climUse(2);

            if climUse(2) == climUse(1)
                colourIdx = ones(size(anglePlot));
            else
                colourIdx = 1 + round( ...
                    (anglePlot - climUse(1)) ./ ...
                    (climUse(2) - climUse(1)) .* (nColours - 1));
            end

            colourIdx = max(1, min(nColours, colourIdx));

            for ic = 1:nColours

                useThisColour = colourIdx == ic;

                if ~any(useThisColour)
                    continue
                end

                idxThis = idxGoodPlot(useThisColour);

                xLines = [xA(idxThis), xB(idxThis), nan(numel(idxThis), 1)]';
                yLines = [yA(idxThis), yB(idxThis), nan(numel(idxThis), 1)]';

                plot(xLines(:), yLines(:), '-', ...
                    'Color', cmap(ic,:), ...
                    'LineWidth', vectorLineWidth);
            end

            colormap(gca, cmap);
            caxis(climUse);

            cb = colorbar;

            if pixelAngleCorrectionDeg == 0
                ylabel(cb, 'Angular residual in pixel-space [deg]');
            else
                ylabel(cb, sprintf('Pixel-space angular deviation - %.1f deg', ...
                    pixelAngleCorrectionDeg));
            end
        end

    else

        xLines = [xA(idxVectorsToPlot), xB(idxVectorsToPlot), nan(numel(idxVectorsToPlot),1)]';
        yLines = [yA(idxVectorsToPlot), yB(idxVectorsToPlot), nan(numel(idxVectorsToPlot),1)]';

        plot(xLines(:), yLines(:), '-', ...
            'Color', [0 1 1], ...
            'LineWidth', vectorLineWidth);
    end
end

if ~isempty(polyUVPlotAll)
    plot(polyUVPlotAll(:,1), polyUVPlotAll(:,2), '-', ...
        'Color', polygonColor, ...
        'LineWidth', polygonLineWidth);
end

% Apply font formatting
set(findall(figH, '-property', 'FontName'), 'FontName', plotFontName);
set(findall(figH, '-property', 'FontSize'), 'FontSize', plotFontSize);

% Re-assert figure size in centimetres
set(figH, 'Units', 'centimeters');
set(figH, 'Position', [figureLeftCm, figureBottomCm, figureWidthCm, figureHeightCm]);
set(figH, 'PaperUnits', 'centimeters');
set(figH, 'PaperPosition', [0, 0, figureWidthCm, figureHeightCm]);
set(figH, 'PaperSize', [figureWidthCm, figureHeightCm]);

%% ------------------------------------------------------------------------
% Save outputs
% -------------------------------------------------------------------------

imageBoundaryUV = polyUVPlotAll;
worldBoundaryXYZ = worldBoundaryXYZAll;
rotatedBoundarySN = rotatedBoundarySNAll;
wseBoundaryColRow = wseBoundaryColRowAll;

wseAdjustedMask = dataMask;

flowVectorIdxInsideAdjustedMask = idxVectorsInRegion;
flowVectorIdxPlotted = idxVectorsToPlot;
flowVectorInsideModeUsed = vectorInsideMode;
flowVectorPixelAngleDeviationPlotted = pixelAngleDeviationPlot;

save(outputFile, ...
    'imageBoundaryUV', ...
    'imageBoundaryUVCell', ...
    'imageBoundaryFiniteMaskCell', ...
    'worldBoundaryXYZ', ...
    'worldBoundaryXYZCell', ...
    'rotatedBoundarySN', ...
    'rotatedBoundarySNCell', ...
    'wseBoundaryColRow', ...
    'wseBoundaryColRowCell', ...
    'wseAdjustedMask', ...
    'imageRegionMask', ...
    'zPlane', ...
    'mapIdx', ...
    'aa', ...
    'flowVectorIdxInsideAdjustedMask', ...
    'flowVectorIdxPlotted', ...
    'flowVectorInsideModeUsed', ...
    'flowVectorPixelAngleDeviationPlotted');

fprintf('\nSaved projected adjusted WSE region and vector selection to:\n  %s\n', outputFile);

exportgraphics(figH, figureFile, 'Resolution', 1200);

outputs = struct( ...
    'figureFile', figureFile, ...
    'diagnosticFile', outputFile, ...
    'inputFile', inputFile, ...
    'checkpointFile', checkpointFile);
end

%% =========================================================================
% Local helper functions
% =========================================================================

function uv = local_get_wse_vector_cell(C, aa, vectorCellIdx, nameText)

    if iscell(C)
        if ndims(C) >= 2 && size(C,1) >= aa && size(C,2) >= vectorCellIdx && ...
                ~isempty(C{aa, vectorCellIdx})

            uv = C{aa, vectorCellIdx};

        elseif numel(C) >= aa && ~isempty(C{aa})

            uv = C{aa};

        else

            error('Could not find %s{%d,%d} or %s{%d}.', ...
                nameText, aa, vectorCellIdx, nameText, aa);
        end
    else
        uv = C;
    end

    if isempty(uv) || size(uv,2) < 2
        error('%s does not contain valid [x y] vector points.', nameText);
    end
end


function inMask = local_image_points_inside_mask(x, y, finitePts, imageMask)

    inMask = false(size(x));

    imageRows = size(imageMask, 1);
    imageCols = size(imageMask, 2);

    good = finitePts & isfinite(x) & isfinite(y);

    col = round(x(good));
    row = round(y(good));

    insideImage = row >= 1 & row <= imageRows & ...
                  col >= 1 & col <= imageCols;

    idxGood = find(good);
    idxInside = idxGood(insideImage);

    rowInside = row(insideImage);
    colInside = col(insideImage);

    lin = sub2ind(size(imageMask), rowInside, colInside);

    inMask(idxInside) = imageMask(lin);
end


function ang = local_wrap180(ang)

    ang = mod(ang + 180, 360) - 180;
end
