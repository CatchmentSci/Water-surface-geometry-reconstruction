function out = make_velocity_maps_from_checkpoint(checkpointInput, inputMatFile, refinedResult, mapTable, opts)
%KLT_MAKE_VELOCITY_MAPS_FROM_CHECKPOINT
% Create flow-velocity maps on the same grid as checkpoint.wse_map.
%
% The function projects KLT image-coordinate path endpoints xyzA_wse/xyzB_wse
% onto two height maps:
%   1) the starting WSE map, normally checkpoint.wse_map{aa,1};
%   2) the accepted/selected WSE map.
%
% It then converts each projected path into velocity in metres per second and
% bins the velocity to the same Ny-by-Nx grid as the WSE map.
%
% Required inputs
%   checkpointInput : checkpoint struct or MAT file containing variable checkpoint
%   inputMatFile    : original solver-input MAT file containing at least:
%                       app_in, xyzA_wse, xyzB_wse
%                     and preferably camA_fullmodel so app.camA can be rebuilt.
%
% Optional inputs
%   refinedResult   : result from KLT_select_wave_checkpoint_refined_first10m
%                     or KLT_select_wave_checkpoint_refined. Used to find the
%                     accepted map if opts.acceptedMapIndex is not supplied.
%   mapTable        : selector map table. Used only if selectedMapIndex is a
%                     row number rather than a checkpoint wse_map index.
%   opts            : struct of options.
%
% Key opts
%   opts.acceptedMapIndex    [] or scalar checkpoint.wse_map index.
%   opts.startMapIndex       default 1.
%   opts.pathPopulation      'opt_idx_full' (default) or 'all'.
%   opts.binAnchor           'A' (default), 'midpoint', or 'B'.
%   opts.dt_s                [] = app.iter/app.videoFrameRate.
%   opts.useRobustMedian     true = median per grid cell, false = mean.
%   opts.minSpeed_mps        default 0, applied after projection.
%   opts.maxSpeed_mps        default Inf, applied after projection.
%   opts.keepProjectedPaths  false; true stores path-level projected data.
%   opts.saveMatFile         '' or filename to save out.
%
% Output fields
%   out.start.speed_mps, out.start.u_streamwise_mps, out.start.v_crossstream_mps
%   out.accepted.speed_mps, out.accepted.u_streamwise_mps, ...
%   out.delta.speed_mps      accepted - start speed map
%   out.grid                 xi, yi, X_rot, Y_rot, etc.
%
% Example
%   S = load('case_ckpt.mat','checkpoint');
%   R = KLT_select_wave_checkpoint_refined_first10m(S.checkpoint);
%   out = KLT_make_velocity_maps_from_checkpoint( ...
%       S.checkpoint, 'case_solver_inputs.mat', R);
%   imagesc(out.accepted.speed_mps); axis image; colorbar

    if nargin < 1 || isempty(checkpointInput)
        if evalin('base', 'exist(''checkpoint'', ''var'')')
            checkpoint = evalin('base', 'checkpoint');
            checkpointFile = '';
        else
            error('Provide a checkpoint struct/file, or define checkpoint in the base workspace.');
        end
    elseif ischar(checkpointInput) || isstring(checkpointInput)
        checkpointFile = char(checkpointInput);
        Sck = load(checkpointFile, 'checkpoint');
        if ~isfield(Sck, 'checkpoint')
            error('checkpointInput MAT file must contain variable checkpoint.');
        end
        checkpoint = Sck.checkpoint;
    elseif isstruct(checkpointInput)
        checkpoint = checkpointInput;
        checkpointFile = '';
    else
        error('checkpointInput must be a checkpoint struct, MAT filename, or empty.');
    end

    if nargin < 2 || isempty(inputMatFile)
        error(['inputMatFile is required. The checkpoint contains WSE maps and grid geometry, ', ...
               'but the raw KLT path endpoints xyzA_wse/xyzB_wse and camera model are normally ', ...
               'stored in the original solver-input MAT file.']);
    end
    if nargin < 3, refinedResult = []; end
    if nargin < 4, mapTable = []; end
    if nargin < 5 || isempty(opts), opts = struct; end
    opts = set_velocity_defaults(opts);

    if isfield(checkpoint, 'aa') && isfinite(checkpoint.aa)
        aa = checkpoint.aa;
    else
        aa = 1;
    end

    requiredCheckpointFields = {'wse_map','xi','yi','X_rot','Y_rot','origin','e_s','e_n'};
    for k = 1:numel(requiredCheckpointFields)
        if ~isfield(checkpoint, requiredCheckpointFields{k})
            error('checkpoint.%s is required.', requiredCheckpointFields{k});
        end
    end

    [wseMaps, mapIndices] = collect_wse_maps(checkpoint.wse_map, aa);
    if isempty(wseMaps)
        error('No non-empty checkpoint.wse_map entries found for aa=%d.', aa);
    end

    startMapIndex = round(opts.startMapIndex);
    acceptedMapIndex = resolve_accepted_map_index(opts, refinedResult, mapTable, mapIndices);

    Wstart = get_wse_map_by_index(checkpoint.wse_map, aa, startMapIndex);
    Wacc = get_wse_map_by_index(checkpoint.wse_map, aa, acceptedMapIndex);

    [Ny, Nx] = size(Wstart);
    if ~isequal(size(Wacc), [Ny Nx])
        error('Starting and accepted WSE maps have different sizes.');
    end

    xi = checkpoint.xi(:).';
    yi = checkpoint.yi(:);
    if numel(xi) ~= Nx || numel(yi) ~= Ny
        error('checkpoint.xi/yi sizes do not match WSE map dimensions.');
    end

    inputData = load_solver_input_for_velocity(inputMatFile);
    app = inputData.app;
    xyzA_wse = inputData.xyzA_wse;
    xyzB_wse = inputData.xyzB_wse;

    if numel(xyzA_wse) < aa || numel(xyzB_wse) < aa || isempty(xyzA_wse{aa}) || isempty(xyzB_wse{aa})
        error('xyzA_wse{%d} and xyzB_wse{%d} must exist and be non-empty.', aa, aa);
    end

    uvA_all = xyzA_wse{aa};
    uvB_all = xyzB_wse{aa};
    if size(uvA_all,1) ~= size(uvB_all,1)
        error('xyzA_wse{%d} and xyzB_wse{%d} have different numbers of paths.', aa, aa);
    end

    pathIdx = choose_path_population(checkpoint, opts, size(uvA_all,1));
    uvA = uvA_all(pathIdx, :);
    uvB = uvB_all(pathIdx, :);

    dt_s = resolve_dt_seconds(app, opts);

    geom = struct;
    geom.xi = xi;
    geom.yi = yi;
    geom.X_rot = checkpoint.X_rot;
    geom.Y_rot = checkpoint.Y_rot;
    geom.origin = checkpoint.origin(:).';
    geom.e_s = checkpoint.e_s(:).';
    geom.e_n = checkpoint.e_n(:).';
    geom.Nx = Nx;
    geom.Ny = Ny;

    startVelocity = velocity_for_height_map(app, uvA, uvB, pathIdx, Wstart, geom, dt_s, opts);
    acceptedVelocity = velocity_for_height_map(app, uvA, uvB, pathIdx, Wacc, geom, dt_s, opts);

    out = struct;
    out.start = startVelocity;
    out.accepted = acceptedVelocity;
    out.delta = struct;
    out.delta.speed_mps = acceptedVelocity.speed_mps - startVelocity.speed_mps;
    out.delta.u_streamwise_mps = acceptedVelocity.u_streamwise_mps - startVelocity.u_streamwise_mps;
    out.delta.v_crossstream_mps = acceptedVelocity.v_crossstream_mps - startVelocity.v_crossstream_mps;
    out.delta.u_world_mps = acceptedVelocity.u_world_mps - startVelocity.u_world_mps;
    out.delta.v_world_mps = acceptedVelocity.v_world_mps - startVelocity.v_world_mps;

    out.grid = geom;
    out.grid.dx_m = median(abs(diff(xi)), 'omitnan');
    out.grid.dy_m = median(abs(diff(yi)), 'omitnan');
    out.aa = aa;
    out.checkpointFile = checkpointFile;
    out.inputMatFile = char(inputMatFile);
    out.startMapIndex = startMapIndex;
    out.acceptedMapIndex = acceptedMapIndex;
    out.dt_s = dt_s;
    out.pathPopulation = opts.pathPopulation;
    out.nInputPaths = size(uvA_all,1);
    out.nUsedPaths = numel(pathIdx);
    out.pathIdx = pathIdx;
    out.binAnchor = opts.binAnchor;
    out.units = 'm/s';
    out.description = ['Velocity maps produced by projecting KLT endpoints onto ', ...
        'the starting WSE map and the accepted WSE map, then binning path velocities ', ...
        'onto the checkpoint.wse_map grid.'];

    if ~isempty(opts.saveMatFile)
        velocityOut = out; %#ok<NASGU>
        save(opts.saveMatFile, 'velocityOut', '-v7.3');
        fprintf('Saved velocity maps to %s\n', opts.saveMatFile);
    end
end

%% ------------------------------------------------------------------------
% Core velocity calculation
%% ------------------------------------------------------------------------
function V = velocity_for_height_map(app, uvA, uvB, pathIdx, H, geom, dt_s, opts)
    H = double(H);

    xyzA = app.camA.invproject(uvA, geom.X_rot, geom.Y_rot, H);
    xyzB = app.camA.invproject(uvB, geom.X_rot, geom.Y_rot, H);

    finite = all(isfinite(xyzA(:,1:2)),2) & all(isfinite(xyzB(:,1:2)),2);

    dxy = xyzB(:,1:2) - xyzA(:,1:2);
    velXY = dxy ./ dt_s;
    speed = hypot(velXY(:,1), velXY(:,2));
    u_s = velXY(:,1).*geom.e_s(1) + velXY(:,2).*geom.e_s(2);
    v_n = velXY(:,1).*geom.e_n(1) + velXY(:,2).*geom.e_n(2);

    finite = finite & isfinite(speed) & speed >= opts.minSpeed_mps & speed <= opts.maxSpeed_mps;

    switch lower(opts.binAnchor)
        case 'a'
            anchorXY = xyzA(:,1:2);
        case 'b'
            anchorXY = xyzB(:,1:2);
        case 'midpoint'
            anchorXY = 0.5 .* (xyzA(:,1:2) + xyzB(:,1:2));
        otherwise
            error('opts.binAnchor must be ''A'', ''B'', or ''midpoint''.');
    end

    [row, col, inGrid] = xy_to_grid_cells(anchorXY(:,1), anchorXY(:,2), geom);
    good = finite & inGrid;
    lin = sub2ind([geom.Ny geom.Nx], row(good), col(good));

    V = struct;
    V.heightMap = H;
    V.speed_mps = bin_values(lin, speed(good), [geom.Ny geom.Nx], opts.useRobustMedian);
    V.u_streamwise_mps = bin_values(lin, u_s(good), [geom.Ny geom.Nx], opts.useRobustMedian);
    V.v_crossstream_mps = bin_values(lin, v_n(good), [geom.Ny geom.Nx], opts.useRobustMedian);
    V.u_world_mps = bin_values(lin, velXY(good,1), [geom.Ny geom.Nx], opts.useRobustMedian);
    V.v_world_mps = bin_values(lin, velXY(good,2), [geom.Ny geom.Nx], opts.useRobustMedian);
    V.count = reshape(accumarray(lin, 1, [geom.Ny*geom.Nx 1], @sum, 0), geom.Ny, geom.Nx);
    V.validPathMask = good;
    V.nProjectedPaths = numel(speed);
    V.nFiniteProjectedPaths = nnz(finite);
    V.nBinnedPaths = nnz(good);
    V.summary = summarize_velocity_map(V.speed_mps, V.count);

    if opts.keepProjectedPaths
        V.paths = struct;
        V.paths.pathIdx = pathIdx(:);
        V.paths.xyzA = xyzA;
        V.paths.xyzB = xyzB;
        V.paths.anchorRow = row;
        V.paths.anchorCol = col;
        V.paths.speed_mps = speed;
        V.paths.u_streamwise_mps = u_s;
        V.paths.v_crossstream_mps = v_n;
        V.paths.u_world_mps = velXY(:,1);
        V.paths.v_world_mps = velXY(:,2);
        V.paths.validPathMask = good;
    end
end

function M = bin_values(lin, vals, gridSize, useMedian)
    nCells = prod(gridSize);
    vals = vals(:);
    good = isfinite(vals) & isfinite(lin) & lin >= 1 & lin <= nCells;
    lin = lin(good);
    vals = vals(good);
    if isempty(lin)
        M = NaN(gridSize);
        return
    end
    if useMedian
        x = accumarray(lin(:), vals(:), [nCells 1], @(v) median(v, 'omitnan'), NaN);
    else
        x = accumarray(lin(:), vals(:), [nCells 1], @(v) mean(v, 'omitnan'), NaN);
    end
    M = reshape(x, gridSize);
end

function [row, col, inGrid] = xy_to_grid_cells(x, y, geom)
    x = x(:);
    y = y(:);
    relX = x - geom.origin(1);
    relY = y - geom.origin(2);
    s = relX .* geom.e_s(1) + relY .* geom.e_s(2);
    n = relX .* geom.e_n(1) + relY .* geom.e_n(2);

    dx = median(abs(diff(geom.xi)), 'omitnan');
    dy = median(abs(diff(geom.yi)), 'omitnan');
    if ~isfinite(dx) || dx <= 0, dx = 0; end
    if ~isfinite(dy) || dy <= 0, dy = 0; end

    col = round(interp1(geom.xi, 1:geom.Nx, s, 'linear', NaN));
    row = round(interp1(geom.yi, 1:geom.Ny, n, 'linear', NaN));

    sMin = min(geom.xi) - 0.5*dx;
    sMax = max(geom.xi) + 0.5*dx;
    nMin = min(geom.yi) - 0.5*dy;
    nMax = max(geom.yi) + 0.5*dy;

    inGrid = isfinite(row) & isfinite(col) & row >= 1 & row <= geom.Ny & col >= 1 & col <= geom.Nx & ...
             s >= sMin & s <= sMax & n >= nMin & n <= nMax;

    row(~inGrid) = 1;
    col(~inGrid) = 1;
    row = double(row);
    col = double(col);
end

%% ------------------------------------------------------------------------
% Input / map helpers
%% ------------------------------------------------------------------------
function opts = set_velocity_defaults(opts)
    opts = default_field(opts, 'startMapIndex', 1);
    opts = default_field(opts, 'acceptedMapIndex', []);
    opts = default_field(opts, 'pathPopulation', 'opt_idx_full');
    opts = default_field(opts, 'binAnchor', 'A');
    opts = default_field(opts, 'dt_s', []);
    opts = default_field(opts, 'useRobustMedian', true);
    opts = default_field(opts, 'minSpeed_mps', 0);
    opts = default_field(opts, 'maxSpeed_mps', Inf);
    opts = default_field(opts, 'keepProjectedPaths', false);
    opts = default_field(opts, 'saveMatFile', '');
end

function s = default_field(s, f, v)
    if ~isfield(s, f) || isempty(s.(f))
        s.(f) = v;
    end
end

function inputData = load_solver_input_for_velocity(inputMatFile)
    if ~(ischar(inputMatFile) || isstring(inputMatFile))
        error('inputMatFile must be a filename.');
    end
    inputMatFile = char(inputMatFile);
    if ~isfile(inputMatFile)
        error('Input MAT file not found: %s', inputMatFile);
    end

    S = load(inputMatFile);
    if ~isfield(S, 'xyzA_wse') || ~isfield(S, 'xyzB_wse')
        error('inputMatFile must contain xyzA_wse and xyzB_wse.');
    end

    if isfield(S, 'app_in')
        app = S.app_in;
    elseif isfield(S, 'app')
        app = S.app;
    else
        error('inputMatFile must contain app_in or app.');
    end

    if isfield(S, 'camA_fullmodel')
        try
            app.camA = camera(S.camA_fullmodel);
        catch ME
            if ~isfield(app, 'camA')
                error('Could not rebuild app.camA from camA_fullmodel: %s', ME.message);
            else
                warning('Could not rebuild app.camA from camA_fullmodel; using app.camA already present. %s');
            end
        end
    elseif ~isfield(app, 'camA')
        error('Could not find camA_fullmodel and app.camA is not present.');
    end

    if isfield(S, 'camA_first_fullmodel')
        try
            app.camA_first = camera(S.camA_first_fullmodel);
        catch
            % camA_first is not needed for velocity projection here.
        end
    end

    inputData = struct;
    inputData.app = app;
    inputData.xyzA_wse = S.xyzA_wse;
    inputData.xyzB_wse = S.xyzB_wse;
end

function pathIdx = choose_path_population(checkpoint, opts, nAll)
    switch lower(opts.pathPopulation)
        case 'opt_idx_full'
            if isfield(checkpoint, 'opt_idx_full') && ~isempty(checkpoint.opt_idx_full)
                pathIdx = double(checkpoint.opt_idx_full(:));
            elseif isfield(checkpoint, 'metrics_static') && isfield(checkpoint.metrics_static, 'optIdxFull')
                pathIdx = double(checkpoint.metrics_static.optIdxFull(:));
            else
                warning('checkpoint.opt_idx_full was not found; using all KLT paths.');
                pathIdx = (1:nAll).';
            end
        case 'all'
            pathIdx = (1:nAll).';
        otherwise
            error('opts.pathPopulation must be ''opt_idx_full'' or ''all''.');
    end
    pathIdx = pathIdx(isfinite(pathIdx) & pathIdx >= 1 & pathIdx <= nAll);
    pathIdx = unique(round(pathIdx), 'stable');
    if isempty(pathIdx)
        error('No valid path indices selected.');
    end
end

function dt_s = resolve_dt_seconds(app, opts)
    if ~isempty(opts.dt_s) && isnumeric(opts.dt_s) && isscalar(opts.dt_s) && isfinite(opts.dt_s) && opts.dt_s > 0
        dt_s = double(opts.dt_s);
        return
    end
    if isfield(app, 'iter') && isfield(app, 'videoFrameRate') && isfinite(app.iter) && isfinite(app.videoFrameRate) && app.videoFrameRate > 0
        dt_s = double(app.iter) ./ double(app.videoFrameRate);
    else
        error('Could not infer dt_s. Set opts.dt_s, or ensure app.iter and app.videoFrameRate are present.');
    end
end

function acceptedMapIndex = resolve_accepted_map_index(opts, refinedResult, mapTable, mapIndices)
    if ~isempty(opts.acceptedMapIndex) && isnumeric(opts.acceptedMapIndex) && isscalar(opts.acceptedMapIndex) && isfinite(opts.acceptedMapIndex)
        acceptedMapIndex = round(opts.acceptedMapIndex);
        return
    end

    acceptedMapIndex = NaN;
    if isstruct(refinedResult) && isfield(refinedResult, 'selectedMapIndex') && isfinite(refinedResult.selectedMapIndex)
        selected = round(refinedResult.selectedMapIndex);
        if ~isempty(mapTable) && istable(mapTable) && ismember('mapIndex', mapTable.Properties.VariableNames) && ...
                selected >= 1 && selected <= height(mapTable)
            acceptedMapIndex = round(mapTable.mapIndex(selected));
        else
            acceptedMapIndex = selected;
        end
    end

    if ~isfinite(acceptedMapIndex)
        acceptedMapIndex = max(mapIndices);
        warning('No acceptedMapIndex/refinedResult supplied; using final non-empty WSE map %d.', acceptedMapIndex);
    end
end

function [maps, mapIndices] = collect_wse_maps(wse_map, aa)
    maps = {};
    mapIndices = [];
    if ~iscell(wse_map)
        error('checkpoint.wse_map must be a cell array.');
    end
    for kk = 1:size(wse_map, 2)
        if size(wse_map,1) >= aa
            W = wse_map{aa, kk};
            if ~isempty(W) && isnumeric(W) && ismatrix(W)
                maps{end+1,1} = double(W); %#ok<AGROW>
                mapIndices(end+1,1) = kk; %#ok<AGROW>
            end
        end
    end
end

function W = get_wse_map_by_index(wse_map, aa, mapIndex)
    if mapIndex < 1 || size(wse_map,1) < aa || size(wse_map,2) < mapIndex || isempty(wse_map{aa,mapIndex})
        error('checkpoint.wse_map{%d,%d} is not available.', aa, mapIndex);
    end
    W = double(wse_map{aa,mapIndex});
end

function summary = summarize_velocity_map(speedMap, countMap)
    v = speedMap(isfinite(speedMap) & countMap > 0);
    summary = struct;
    summary.nCells = numel(v);
    if isempty(v)
        summary.medianSpeed_mps = NaN;
        summary.meanSpeed_mps = NaN;
        summary.q25Speed_mps = NaN;
        summary.q75Speed_mps = NaN;
        summary.maxSpeed_mps = NaN;
    else
        summary.medianSpeed_mps = median(v, 'omitnan');
        summary.meanSpeed_mps = mean(v, 'omitnan');
        summary.q25Speed_mps = prctile(v, 25);
        summary.q75Speed_mps = prctile(v, 75);
        summary.maxSpeed_mps = max(v);
    end
end
