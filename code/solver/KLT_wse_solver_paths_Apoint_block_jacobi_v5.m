function [wse_map, metrics, phi] = KLT_wse_solver_paths_Apoint_block_jacobi_v5( ...
    app, xyzA_wse, xyzB_wse, aa, wse_map, globalPolarity, ...
    checkpointFile, resumeFromCheckpoint, inputMatFile, ...
    sweepConfigFile, videoNumberForSweep, flowAngleOffsetDeg, ...
    residualZeroRowRange) %#ok<INUSD>
%KLT_WSE_SOLVER_PATHS_APOINT_BLOCK_JACOBI_V3_CONVERGENCE
%
% Column-wise Jacobi/Gauss-Seidel swept-cell WSE solver with A-point
% representative assignment and extended convergence diagnostics.
%
% This version keeps the WSE update logic unchanged relative to the supplied
% v3 solver, except for additional bookkeeping needed to evaluate:
%   1) within-sweep local residual improvement;
%   2) true end-of-iteration residuals from a full reprojection;
%   3) fixed-population path residuals;
%   4) fixed representative-path residuals;
%   5) WSE surface update magnitudes;
%   6) active-set churn;
%   7) update sign oscillation;
%   8) rejection / acceptance reasons for each eligible cell;
%   9) row-band residual and update behaviour;
%  10) diagnostic convergence/stall flags.
%
% No automatic stopping rule is applied.

    % -------------------------------------------------------------------------
    % Optional inputs
    % -------------------------------------------------------------------------
    if nargin < 9
        inputMatFile = [];
    end
    if nargin < 10
        sweepConfigFile = [];
    end
    if nargin < 11
        videoNumberForSweep = [];
    end
    if nargin < 12
        flowAngleOffsetDeg = [];
    end
    if nargin < 13
        residualZeroRowRange = [];
    end

    % The column range is selected internally. When a sweep config is used,
    % both the configured rows and columns define the residual-angle window.
    residualZeroColRange = [];
    flowAngleOffsetWindowSource = 'sweepConfigFile';
    flowAngleOffsetMode = 'sweep_config_window_inferred';

    flowAngleOffsetDegWasInput = ~isempty(flowAngleOffsetDeg) && ...
        isnumeric(flowAngleOffsetDeg) && isscalar(flowAngleOffsetDeg) && ...
        isfinite(double(flowAngleOffsetDeg));
    if flowAngleOffsetDegWasInput
        flowAngleOffsetDeg = double(flowAngleOffsetDeg);
    else
        flowAngleOffsetDeg = [];
    end

    % -------------------------------------------------------------------------
    % Optional load from MAT file
    % If inputMatFile is supplied, it overrides the first 6 inputs.
    % -------------------------------------------------------------------------
    if ~isempty(inputMatFile)
        if ~isfile(inputMatFile)
            error('Input MAT file not found: %s', inputMatFile);
        end

        S = load(inputMatFile);
        requiredVars = { ...
            'app_in', ...
            'camA_fullmodel', ...
            'camA_first_fullmodel', ...
            'xyzA_wse', ...
            'xyzB_wse', ...
            'aa', ...
            'wse_map'};

        for k = 1:numel(requiredVars)
            if ~isfield(S, requiredVars{k})
                error('Input MAT file is missing required variable: %s', requiredVars{k});
            end
        end

        app = S.app_in;
        app.camA = camera(S.camA_fullmodel);
        app.camA_first = camera(S.camA_first_fullmodel);
        xyzA_wse = S.xyzA_wse;
        xyzB_wse = S.xyzB_wse;
        aa = S.aa;
        wse_map = S.wse_map;

        if isfield(S, 'globalPolarity')
            globalPolarity = S.globalPolarity; %#ok<NASGU>
        else
            globalPolarity = []; %#ok<NASGU>
        end
    end

    % -------------------------------------------------------------------------
    % Defaults for checkpoint / restart
    % -------------------------------------------------------------------------
    if nargin < 7 || isempty(checkpointFile)
        checkpointFile = sprintf('KLT_WSE_Apoint_jacobi_convergence_checkpoint_aa_%d.mat', aa);
    end
    if nargin < 8 || isempty(resumeFromCheckpoint)
        resumeFromCheckpoint = false;
    end

    % -------------------------------------------------------------------------
    % Basic input checks
    % -------------------------------------------------------------------------
    if isempty(app)
        error('app is empty.');
    end
    if isempty(xyzA_wse)
        error('xyzA_wse is empty.');
    end
    if isempty(xyzB_wse)
        error('xyzB_wse is empty.');
    end
    if isempty(aa)
        error('aa is empty.');
    end
    if isempty(wse_map)
        error('wse_map is empty.');
    end

    % -------------------------------------------------------------------------
    % Settings
    % -------------------------------------------------------------------------
    increment = 0.005;      % metres; finite-difference perturbation
    alpha = 0.5;           % damping on accepted solved update
    maxIter = 200;
    useParallelRows = KLT_can_use_parallel();
    verbosePolarityPrint = false;

    % Column-wise Jacobi/Newton controls
    jacobiMaxStep = 1 * increment;        % was 0.05 metres; max undamped dz_required works better at  2 x increment with 0.1 being increment
    jacobiMinSensitivity = 1e-6;          % deg/m
    jacobiRequireImprovement = true;
    jacobiFallbackToIncrement = false;

    % Reprojection acceleration switches
    usePreparedDEMSeed = true;
    useWarmStarts = false;
    usePreparedDEMRefresh = true;
    useWarmStartsRefresh = false;
    usePreparedDEMRowTrials = true;

    rebuildMembershipAfterChangedColumn = true; % was false for the benchmark runs
    debugDiagnostics = true;
    enforceCurrentDiffCol = true;
    useStencilAffectedRefresh = false;
    nNearDebugRows = 20;

    % Convergence diagnostic controls
    doEndIterationFullReprojection = true;
    nConvergenceBands = 5;
    convergenceWindow = 5;

    polarityMode = 'jacobi_column';
    metrics = struct;

    % =========================================================================
    % Resume from checkpoint OR do a fresh initialisation
    % =========================================================================
    if resumeFromCheckpoint
        if ~isfile(checkpointFile)
            error('Requested resumeFromCheckpoint=true, but checkpoint file not found: %s', checkpointFile);
        end

        S = load(checkpointFile, 'checkpoint');
        checkpoint = S.checkpoint;

        if checkpoint.aa ~= aa
            error('Checkpoint aa (%d) does not match requested aa (%d).', ...
                checkpoint.aa, aa);
        end

        fprintf('\n============================================================\n');
        fprintf('Resuming from checkpoint: %s\n', checkpointFile);
        fprintf('Last completed iteration: %d\n', checkpoint.lastCompletedIter);
        fprintf('Next iteration to run: %d\n', checkpoint.nextIter);
        fprintf('============================================================\n');

        % Restore solver settings
        increment = checkpoint.increment;
        alpha = checkpoint.alpha;
        if maxIter < checkpoint.maxIter
            maxIter = checkpoint.maxIter;
        end
        polarityMode = checkpoint.polarityMode;
        if isfield(checkpoint, 'verbosePolarityPrint')
            verbosePolarityPrint = checkpoint.verbosePolarityPrint;
        end

        if isfield(checkpoint, 'jacobiMaxStep')
            jacobiMaxStep = checkpoint.jacobiMaxStep;
        else
            jacobiMaxStep = 5 * increment;
        end
        if isfield(checkpoint, 'jacobiMinSensitivity')
            jacobiMinSensitivity = checkpoint.jacobiMinSensitivity;
        else
            jacobiMinSensitivity = 1e-6;
        end
        if isfield(checkpoint, 'jacobiRequireImprovement')
            jacobiRequireImprovement = checkpoint.jacobiRequireImprovement;
        else
            jacobiRequireImprovement = true;
        end
        if isfield(checkpoint, 'jacobiFallbackToIncrement')
            jacobiFallbackToIncrement = checkpoint.jacobiFallbackToIncrement;
        else
            jacobiFallbackToIncrement = false;
        end

        % Restore geometry / static state
        phi = checkpoint.phi;
        shift = checkpoint.shift;
        e_s = checkpoint.e_s;
        e_n = checkpoint.e_n;
        origin = checkpoint.origin;
        xi = checkpoint.xi;
        yi = checkpoint.yi;
        X_rot = checkpoint.X_rot;
        Y_rot = checkpoint.Y_rot;
        Nx = checkpoint.Nx;
        Ny = checkpoint.Ny;
        dx = checkpoint.dx;
        dy = checkpoint.dy;
        x_order = checkpoint.x_order;
        y_order = checkpoint.y_order;
        active_window_mask = checkpoint.active_window_mask;
        nSweepCells = checkpoint.nSweepCells;
        opt_idx_full = checkpoint.opt_idx_full;
        cell_selected_idx_init = checkpoint.cell_selected_idx_init;
        cell_median_val_init = checkpoint.cell_median_val_init;

        % Restore dynamic state
        wse_map = checkpoint.wse_map;
        idx = checkpoint.idx;
        startIter = checkpoint.nextIter;

        % Restore diagnostics
        zi = KLT_checkpoint_cell_or_default(checkpoint, 'zi', maxIter);
        zi_post = KLT_checkpoint_cell_or_default(checkpoint, 'zi_post', maxIter);
        dzMapHist = KLT_checkpoint_cell_or_default(checkpoint, 'dzMapHist', maxIter);
        signMapHist = KLT_checkpoint_cell_or_default(checkpoint, 'signMapHist', maxIter);
        residualMapHist = KLT_checkpoint_cell_or_default(checkpoint, 'residualMapHist', maxIter);
        residualPostMapHist = KLT_checkpoint_cell_or_default(checkpoint, 'residualPostMapHist', maxIter);
        selectedPathMapHist = KLT_checkpoint_cell_or_default(checkpoint, 'selectedPathMapHist', maxIter);
        xi_hist = KLT_checkpoint_cell_or_default(checkpoint, 'xi_hist', maxIter);
        yi_hist = KLT_checkpoint_cell_or_default(checkpoint, 'yi_hist', maxIter);

        meanAbsResidualHist = KLT_checkpoint_vec_or_default(checkpoint, 'meanAbsResidualHist', maxIter);
        medianAbsResidualHist = KLT_checkpoint_vec_or_default(checkpoint, 'medianAbsResidualHist', maxIter);
        rmsResidualHist = KLT_checkpoint_vec_or_default(checkpoint, 'rmsResidualHist', maxIter);
        meanSignedResidualHist = KLT_checkpoint_vec_or_default(checkpoint, 'meanSignedResidualHist', maxIter);
        medianSignedResidualHist = KLT_checkpoint_vec_or_default(checkpoint, 'medianSignedResidualHist', maxIter);

        meanAbsResidualPostHist = KLT_checkpoint_vec_or_default(checkpoint, 'meanAbsResidualPostHist', maxIter);
        medianAbsResidualPostHist = KLT_checkpoint_vec_or_default(checkpoint, 'medianAbsResidualPostHist', maxIter);
        rmsResidualPostHist = KLT_checkpoint_vec_or_default(checkpoint, 'rmsResidualPostHist', maxIter);
        meanSignedResidualPostHist = KLT_checkpoint_vec_or_default(checkpoint, 'meanSignedResidualPostHist', maxIter);
        medianSignedResidualPostHist = KLT_checkpoint_vec_or_default(checkpoint, 'medianSignedResidualPostHist', maxIter);

        medianAbsDzHist = KLT_checkpoint_vec_or_default(checkpoint, 'medianAbsDzHist', maxIter);
        maxAbsDzHist = KLT_checkpoint_vec_or_default(checkpoint, 'maxAbsDzHist', maxIter);

        activeCellFracHist = KLT_checkpoint_vec_or_default(checkpoint, 'activeCellFracHist', maxIter);
        updatedCellFracHist = KLT_checkpoint_vec_or_default(checkpoint, 'updatedCellFracHist', maxIter);

        if isfield(checkpoint, 'metrics_static')
            metrics = checkpoint.metrics_static;
        else
            metrics = struct;
        end

        if isempty(residualZeroRowRange)
            if isfield(checkpoint, 'residualZeroRowRange')
                residualZeroRowRange = checkpoint.residualZeroRowRange;
            elseif isfield(metrics, 'residualZeroRowRange')
                residualZeroRowRange = metrics.residualZeroRowRange;
            end
        end

        if isfield(checkpoint, 'residualZeroColRange')
            residualZeroColRange = checkpoint.residualZeroColRange;
        elseif isfield(metrics, 'residualZeroColRange')
            residualZeroColRange = metrics.residualZeroColRange;
        end

        % Preserve the sweep-window residual-angle calibration on resume.
        % A newly supplied config or a config recorded in the checkpoint means
        % that the restored sweep rows and columns define the calibration area.
        useSweepWindowForFlowAngleOffset = ~isempty(sweepConfigFile);
        if ~useSweepWindowForFlowAngleOffset && ...
                isfield(metrics, 'sweepConfigFile') && ~isempty(metrics.sweepConfigFile)
            useSweepWindowForFlowAngleOffset = true;
        end

        if useSweepWindowForFlowAngleOffset
            residualZeroRowRange = y_order;
            residualZeroColRange = x_order;
            flowAngleOffsetWindowSource = 'sweepConfigFile';
        elseif ~isempty(residualZeroRowRange)
            residualZeroColRange = [];
            flowAngleOffsetWindowSource = 'residualZeroRowRange';
        else
            residualZeroColRange = [];
            flowAngleOffsetWindowSource = 'full_domain';
        end

        if isfield(checkpoint, 'residualZeroMedianRawDeg')
            residualZeroMedianRawDeg = checkpoint.residualZeroMedianRawDeg;
        elseif isfield(metrics, 'residualZeroMedianRawDeg')
            residualZeroMedianRawDeg = metrics.residualZeroMedianRawDeg;
        else
            residualZeroMedianRawDeg = NaN;
        end

        if useSweepWindowForFlowAngleOffset
            % A supplied sweep config takes precedence over a manually passed
            % angle offset: infer the offset from its configured rows/columns.
            [flowAngleOffsetDeg, residualZeroMedianRawDeg] = ...
                KLT_infer_flow_angle_offset_from_metrics(metrics, ...
                residualZeroRowRange, residualZeroColRange);
            flowAngleOffsetMode = 'sweep_config_window_inferred';
        elseif flowAngleOffsetDegWasInput
            flowAngleOffsetDeg = double(flowAngleOffsetDeg);
            flowAngleOffsetMode = 'explicit_input';
            if isfield(metrics, 'initCellMedianValRaw') && ~isempty(metrics.initCellMedianValRaw)
                [~, residualZeroMedianRawDeg] = ...
                    KLT_median_residual_in_window(metrics.initCellMedianValRaw, ...
                    residualZeroRowRange, residualZeroColRange);
            elseif isfield(metrics, 'initCellMedianVal') && ~isempty(metrics.initCellMedianVal)
                [~, residualZeroMedianRawDeg] = ...
                    KLT_median_residual_in_window(metrics.initCellMedianVal, ...
                    residualZeroRowRange, residualZeroColRange);
            end
        else
            [flowAngleOffsetDeg, residualZeroMedianRawDeg] = ...
                KLT_infer_flow_angle_offset_from_metrics(metrics, ...
                residualZeroRowRange, residualZeroColRange);
            flowAngleOffsetMode = 'residual_median_inferred';
        end

        if isempty(flowAngleOffsetDeg) || ~isfinite(flowAngleOffsetDeg)
            warning('Residual-angle offset was not finite during checkpoint resume; using 0 degrees.');
            flowAngleOffsetDeg = 0;
        end

        metrics.flowAngleOffsetDeg = flowAngleOffsetDeg;
        metrics.flowAngleOffsetDegWasInput = flowAngleOffsetDegWasInput;
        metrics.flowAngleOffsetWindowSource = flowAngleOffsetWindowSource;
        metrics.flowAngleOffsetMode = flowAngleOffsetMode;
        metrics.residualZeroRowRange = residualZeroRowRange;
        metrics.residualZeroColRange = residualZeroColRange;
        metrics.residualZeroMedianRawDeg = residualZeroMedianRawDeg;
        metrics.jacobiMaxStep = jacobiMaxStep;
        metrics.jacobiMinSensitivity = jacobiMinSensitivity;
        metrics.jacobiRequireImprovement = jacobiRequireImprovement;
        metrics.jacobiFallbackToIncrement = jacobiFallbackToIncrement;

        app.Transdem = wse_map{aa, idx};

    else
        % ---------------------------------------------------------------------
        % Diagnostics storage
        % ---------------------------------------------------------------------
        zi = cell(maxIter,1);
        zi_post = cell(maxIter,1);
        dzMapHist = cell(maxIter,1);
        signMapHist = cell(maxIter,1);
        residualMapHist = cell(maxIter,1);
        residualPostMapHist = cell(maxIter,1);
        selectedPathMapHist = cell(maxIter,1);
        xi_hist = cell(maxIter,1);
        yi_hist = cell(maxIter,1);

        meanAbsResidualHist = nan(maxIter,1);
        medianAbsResidualHist = nan(maxIter,1);
        rmsResidualHist = nan(maxIter,1);
        meanSignedResidualHist = nan(maxIter,1);
        medianSignedResidualHist = nan(maxIter,1);

        meanAbsResidualPostHist = nan(maxIter,1);
        medianAbsResidualPostHist = nan(maxIter,1);
        rmsResidualPostHist = nan(maxIter,1);
        meanSignedResidualPostHist = nan(maxIter,1);
        medianSignedResidualPostHist = nan(maxIter,1);

        medianAbsDzHist = nan(maxIter,1);
        maxAbsDzHist = nan(maxIter,1);

        activeCellFracHist = nan(maxIter,1);
        updatedCellFracHist = nan(maxIter,1);

        % ---------------------------------------------------------------------
        % Ideal downstream direction
        % ---------------------------------------------------------------------
        xIn = app.pts(1,:)';
        yIn = app.pts(2,:)';
        app.start1 = [xIn(1), yIn(1)];
        app.end1 = [xIn(2), yIn(2)];

        frameSize = size(app.firstFrame);
        t1 = frameSize(1);
        t2 = frameSize(2);

        app.start1(1) = max(0, min(t2, app.start1(1)));
        app.start1(2) = max(0, min(t1, app.start1(2)));
        app.end1(1) = max(0, min(t2, app.end1(1)));
        app.end1(2) = max(0, min(t1, app.end1(2)));

        hgt_use = app.riverLevelAnalysis(app.videoNumber);

        if ~strcmp(app.OrientationDropDown.Value, 'Dynamic: Stabilisation') && ...
                ~strcmp(app.OrientationDropDown.Value, 'Planet [beta]')
            tempDEM = zeros(size(app.TransX)) + hgt_use;
            start1_rw_raw = app.camA_first.invproject( ...
                app.start1, app.TransX, app.TransY, tempDEM);
            end1_rw_raw = app.camA_first.invproject( ...
                app.end1, app.TransX, app.TransY, tempDEM);
            start1_rw = KLT_force_xy_point(start1_rw_raw);
            end1_rw = KLT_force_xy_point(end1_rw_raw);
        else
            start1_rw = app.start1 .* app.imageResolution;
            end1_rw = app.end1 .* app.imageResolution;
            start1_rw = start1_rw(1:2);
            end1_rw = end1_rw(1:2);
        end

        ideal = end1_rw - start1_rw;
        ideal = ideal(:).';
        L_ideal = hypot(ideal(1), ideal(2));
        if L_ideal == 0 || ~isfinite(L_ideal)
            error('Ideal downstream direction is invalid.');
        end

        e_s = ideal / L_ideal;
        e_n = [-e_s(2), e_s(1)];

        [ang_out, shift] = KLT_wrapTo360_centerMedian( ...
            rad2deg(atan2(ideal(2), ideal(1))));
        phi = ang_out;

        % ---------------------------------------------------------------------
        % Grid geometry: structured in (s,n), mapped back to world as X_rot/Y_rot
        % ---------------------------------------------------------------------
        dx_native = median(abs(diff(app.X(1,:))), 'omitnan');
        dy_native = median(abs(diff(app.Y(:,1))), 'omitnan');
        target_dx = 0.1;
        target_dy = 0.1;

        allXY = [app.X(:), app.Y(:)];
        good = all(isfinite(allXY), 2);
        allXY = allXY(good, :);

        origin = start1_rw(:).';
        relXY = allXY - origin;
        s_all = relXY * e_s(:);
        n_all = relXY * e_n(:);

        xi = floor(min(s_all)/target_dx)*target_dx : ...
            target_dx : ...
            ceil(max(s_all)/target_dx)*target_dx;
        yi = floor(min(n_all)/target_dy)*target_dy : ...
            target_dy : ...
            ceil(max(n_all)/target_dy)*target_dy;

        [Sg, Ng] = meshgrid(xi, yi);
        X_rot = origin(1) + Sg*e_s(1) + Ng*e_n(1);
        Y_rot = origin(2) + Sg*e_s(2) + Ng*e_n(2);
        Nx = numel(xi);
        Ny = numel(yi);

        % ---------------------------------------------------------------------
        % Sweep region selection
        % ---------------------------------------------------------------------
        analyse_all = false;
        if analyse_all
            col_start = 1;
            col_end = Nx;
            row_near = Ny;
            row_far = 1;
        else
            if ~isempty(sweepConfigFile)
                if isempty(videoNumberForSweep)
                    error(['sweepConfigFile was provided, but videoNumberForSweep is empty. ', ...
                        'Pass the video number explicitly when calling the solver.']);
                end
                [col_start, col_end, row_near, row_far] = ...
                    KLT_load_sweep_limits_from_text(sweepConfigFile, videoNumberForSweep);
            else
                col_start = 1;
                col_end = Nx;
                row_near = Ny;
                row_far = 1;
            end
        end

        dx = target_dx;
        dy = target_dy;

        metrics.dx_native_rot = dx_native;
        metrics.dy_native_rot = dy_native;
        metrics.dx_used = dx;
        metrics.dy_used = dy;
        metrics.e_s = e_s;
        metrics.e_n = e_n;
        metrics.grid_origin = origin;
        metrics.X_rot = X_rot;
        metrics.Y_rot = Y_rot;
        metrics.phi = phi;
        metrics.shift = shift;
        metrics.start1_rw = start1_rw;
        metrics.end1_rw = end1_rw;
        metrics.sweepConfigFile = sweepConfigFile;
        metrics.sweepVideoNumber = videoNumberForSweep;
        metrics.jacobiMaxStep = jacobiMaxStep;
        metrics.jacobiMinSensitivity = jacobiMinSensitivity;
        metrics.jacobiRequireImprovement = jacobiRequireImprovement;
        metrics.jacobiFallbackToIncrement = jacobiFallbackToIncrement;

        % ---------------------------------------------------------------------
        % Initialise flat DEM
        % ---------------------------------------------------------------------
        sizer = size(X_rot);
        app.Transdem = zeros(sizer);
        app.Transdem(:,:) = hgt_use;
        wse_map{aa,1} = app.Transdem;
        idx = 1;

        % ---------------------------------------------------------------------
        % Sweep order
        % ---------------------------------------------------------------------
        col_start = max(1, min(Nx, col_start));
        col_end = max(1, min(Nx, col_end));
        row_near = max(1, min(Ny, row_near));
        row_far = max(1, min(Ny, row_far));

        if col_start <= col_end
            x_order = col_start:col_end;
        else
            x_order = col_start:-1:col_end;
        end

        if row_near >= row_far
            y_order = row_near:-1:row_far;
        else
            y_order = row_near:row_far;
        end

        if isempty(x_order)
            error('x_order is empty after applying column limits.');
        end
        if isempty(y_order)
            error('y_order is empty after applying row limits.');
        end

        % If a sweep config is present, use its effective (clamped) sweep
        % rows and columns to infer the residual-angle flow offset. Without a
        % config, retain the legacy optional row-only band or full-domain mode.
        useSweepWindowForFlowAngleOffset = ~isempty(sweepConfigFile);
        if useSweepWindowForFlowAngleOffset
            residualZeroRowRange = y_order;
            residualZeroColRange = x_order;
            flowAngleOffsetWindowSource = 'sweepConfigFile';
        elseif ~isempty(residualZeroRowRange)
            residualZeroColRange = [];
            flowAngleOffsetWindowSource = 'residualZeroRowRange';
        else
            residualZeroColRange = [];
            flowAngleOffsetWindowSource = 'full_domain';
        end

        metrics.sweepCols = x_order;
        metrics.sweepRows = y_order;
        metrics.sweepX = xi(x_order);
        metrics.sweepY = yi(y_order);
        metrics.polarityMode = polarityMode;

        sweepMask = nan(Ny, Nx);
        sweepMask(y_order, x_order) = 1;
        metrics.sweepMask = sweepMask;

        sweepOrderMap = nan(Ny, Nx);
        kk = 1;
        for ix_tmp = x_order
            for iy_tmp = y_order
                sweepOrderMap(iy_tmp, ix_tmp) = kk;
                kk = kk + 1;
            end
        end
        metrics.sweepOrderMap = sweepOrderMap;

        active_window_mask = false(Ny, Nx);
        active_window_mask(y_order, x_order) = true;
        nSweepCells = nnz(active_window_mask);

        % ---------------------------------------------------------------------
        % Build fixed optimisation subset from initial DEM
        % ---------------------------------------------------------------------
        xyzA0_all = app.camA.invproject( ...
            xyzA_wse{aa}, X_rot, Y_rot, wse_map{aa,1});
        xyzB0_all = app.camA.invproject( ...
            xyzB_wse{aa}, X_rot, Y_rot, wse_map{aa,1});

        vel0_all = xyzB0_all - xyzA0_all;
        dvel0_all = vel0_all ./ (app.iter * 1 / app.videoFrameRate);
        vmag0_all = hypot(dvel0_all(:,1), dvel0_all(:,2));
        keep_speed = vmag0_all > 0.1;
        idx_speed = find(keep_speed);

        v20 = xyzB0_all(idx_speed,1:2) - xyzA0_all(idx_speed,1:2);
        obs_dir0_raw = KLT_applyAngleShift(rad2deg(atan2(v20(:,2), v20(:,1))), shift);
        psi0_raw = phi - obs_dir0_raw;

        filterAngle0 = 90;
        keep_angle = psi0_raw >= -filterAngle0 & psi0_raw <= filterAngle0;
        idx_sa = idx_speed(keep_angle);
        psi_sa_raw = psi0_raw(keep_angle);

        [cellA_x0, cellA_y0] = KLT_points_to_cells( ...
            xyzA0_all(idx_sa,1), xyzA0_all(idx_sa,2), ...
            xi, yi, dx, dy, origin, e_s, e_n);
        [cellB_x0, cellB_y0] = KLT_points_to_cells( ...
            xyzB0_all(idx_sa,1), xyzB0_all(idx_sa,2), ...
            xi, yi, dx, dy, origin, e_s, e_n);

        valid_bins0 = isfinite(cellA_x0) & isfinite(cellA_y0) & ...
            isfinite(cellB_x0) & isfinite(cellB_y0);
        keep_diffcol0 = valid_bins0 & (cellA_x0 ~= cellB_x0);

        idx_cand_full = idx_sa(keep_diffcol0);
        r_cand_raw = psi_sa_raw(keep_diffcol0);
        cand_valid_mask = true(numel(r_cand_raw), 1);

        [~, ~, cell_median_val_init_raw] = ...
            KLT_select_one_path_per_cell( ...
            xyzA0_all(idx_cand_full,1), xyzA0_all(idx_cand_full,2), ...
            r_cand_raw, cand_valid_mask, xi, yi, dx, dy, origin, e_s, e_n);

        metrics.initCellMedianVal = cell_median_val_init_raw;

        if useSweepWindowForFlowAngleOffset
            % A supplied sweep config takes precedence over a manually passed
            % angle offset: infer the offset from its configured rows/columns.
            [flowAngleOffsetDeg, residualZeroMedianRawDeg] = ...
                KLT_infer_flow_angle_offset_from_metrics(metrics, ...
                residualZeroRowRange, residualZeroColRange);
            flowAngleOffsetMode = 'sweep_config_window_inferred';
        elseif flowAngleOffsetDegWasInput
            flowAngleOffsetDeg = double(flowAngleOffsetDeg);
            flowAngleOffsetMode = 'explicit_input';
            [~, residualZeroMedianRawDeg] = ...
                KLT_median_residual_in_window(metrics.initCellMedianVal, ...
                residualZeroRowRange, residualZeroColRange);
        else
            [flowAngleOffsetDeg, residualZeroMedianRawDeg] = ...
                KLT_infer_flow_angle_offset_from_metrics(metrics, ...
                residualZeroRowRange, residualZeroColRange);
            flowAngleOffsetMode = 'residual_median_inferred';
        end

        if isempty(flowAngleOffsetDeg) || ~isfinite(flowAngleOffsetDeg)
            warning('Residual-angle offset was not finite; using 0 degrees.');
            flowAngleOffsetDeg = 0;
        end

        r_cand = r_cand_raw - flowAngleOffsetDeg;

        [selected_rep_mask_init, cell_selected_idx_init, cell_median_val_init] = ...
            KLT_select_one_path_per_cell( ...
            xyzA0_all(idx_cand_full,1), xyzA0_all(idx_cand_full,2), ...
            r_cand, cand_valid_mask, xi, yi, dx, dy, origin, e_s, e_n);

        opt_idx_full = idx_cand_full;

        if isempty(opt_idx_full)
            error('No candidate vectors survived into the fixed WSE optimisation subset.');
        end

        metrics.optIdxFull = opt_idx_full;
        metrics.nOptVectors = numel(opt_idx_full);
        metrics.initRepresentativeMaskWithinCandidatePool = selected_rep_mask_init;
        metrics.initRepresentativeIdxFull = idx_cand_full(selected_rep_mask_init);
        metrics.nInitRepresentativeVectors = nnz(selected_rep_mask_init);
        metrics.initCellChosenIdx = cell_selected_idx_init;
        metrics.initCellMedianVal = cell_median_val_init;
        metrics.initCellMedianValRaw = cell_median_val_init_raw;
        metrics.flowAngleOffsetDeg = flowAngleOffsetDeg;
        metrics.flowAngleOffsetDegWasInput = flowAngleOffsetDegWasInput;
        metrics.flowAngleOffsetWindowSource = flowAngleOffsetWindowSource;
        metrics.flowAngleOffsetMode = flowAngleOffsetMode;
        metrics.residualZeroRowRange = residualZeroRowRange;
        metrics.residualZeroColRange = residualZeroColRange;
        metrics.residualZeroMedianRawDeg = residualZeroMedianRawDeg;

        fprintf('Residual-angle flow offset: %+0.6f deg', flowAngleOffsetDeg);
        if ~isempty(residualZeroRowRange) && ~isempty(residualZeroColRange) && ...
                isfinite(residualZeroMedianRawDeg)
            fprintf(' (raw median %+0.6f deg over rows %d:%d, cols %d:%d)', ...
                residualZeroMedianRawDeg, ...
                min(residualZeroRowRange), max(residualZeroRowRange), ...
                min(residualZeroColRange), max(residualZeroColRange));
        elseif ~isempty(residualZeroRowRange) && isfinite(residualZeroMedianRawDeg)
            fprintf(' (raw median %+0.6f deg over rows %d:%d)', ...
                residualZeroMedianRawDeg, min(residualZeroRowRange), max(residualZeroRowRange));
        elseif ~flowAngleOffsetDegWasInput && isfinite(residualZeroMedianRawDeg)
            fprintf(' (raw median %+0.6f deg across full domain)', ...
                residualZeroMedianRawDeg);
        end
        fprintf('\n');

        startIter = 1;
    end

    if startIter > maxIter
        fprintf('Checkpoint already completed all %d iterations. Returning saved state.\n', maxIter);
    end

    xyzA_opt_uv = xyzA_wse{aa}(opt_idx_full, :);
    xyzB_opt_uv = xyzB_wse{aa}(opt_idx_full, :);

    % Prepared-DEM state
    demPrepBase = app.camA.prepareDEMInverse(X_rot, Y_rot);
    demPrepCurrent = app.camA.updatePreparedDEM(demPrepBase, wse_map{aa,idx});

    prev_xyA_opt = [];
    prev_xyB_opt = [];

    % Debug and advanced diagnostic storage
    if debugDiagnostics
        metrics = KLT_ensure_debug_metrics(metrics, maxIter);

        nNearUse = min(nNearDebugRows, numel(y_order));
        nearBand = double(y_order(1:nNearUse));

        metrics.debug.nearBand = nearBand;
        metrics.debug.nNearDebugRows = nNearUse;
        metrics.debug.enforceCurrentDiffCol = enforceCurrentDiffCol;
        metrics.debug.useStencilAffectedRefresh = useStencilAffectedRefresh;
    else
        nearBand = [];
    end

    metrics = KLT_ensure_advanced_convergence_metrics(metrics, maxIter);
    metrics.convergenceSettings.doEndIterationFullReprojection = doEndIterationFullReprojection;
    metrics.convergenceSettings.nConvergenceBands = nConvergenceBands;
    metrics.convergenceSettings.convergenceWindow = convergenceWindow;

    % =========================================================================
    % Outer iterations
    % =========================================================================
    for a = startIter:maxIter
        fprintf('\n============================================================\n');
        fprintf('Iteration %d / %d\n', a, maxIter);
        fprintf('============================================================\n');

        wse_map{aa, idx+1} = wse_map{aa, idx};
        idx = idx + 1;

        % ---------------------------------------------------------------------
        % Step 1: full reprojection at iteration seed
        % ---------------------------------------------------------------------
        if useWarmStarts && ~isempty(prev_xyA_opt) && ~isempty(prev_xyB_opt) && ...
                size(prev_xyA_opt,1) == size(xyzA_opt_uv,1) && ...
                size(prev_xyB_opt,1) == size(xyzB_opt_uv,1)
            [xyzA_cache, xyzB_cache] = KLT_invproject_uv_pairs_safe_or_prepared( ...
                app.camA, ...
                xyzA_opt_uv, xyzB_opt_uv, ...
                X_rot, Y_rot, wse_map{aa,idx}, ...
                usePreparedDEMSeed, demPrepCurrent, ...
                prev_xyA_opt, prev_xyB_opt);
        else
            [xyzA_cache, xyzB_cache] = KLT_invproject_uv_pairs_safe_or_prepared( ...
                app.camA, ...
                xyzA_opt_uv, xyzB_opt_uv, ...
                X_rot, Y_rot, wse_map{aa,idx}, ...
                usePreparedDEMSeed, demPrepCurrent);
        end

        psi_cache = KLT_compute_residual_from_projected( ...
            xyzA_cache(:,1:2), xyzB_cache(:,1:2), phi, shift, flowAngleOffsetDeg);

        if debugDiagnostics
            finiteA_xy = all(isfinite(xyzA_cache(:,1:2)), 2);
            finiteB_xy = all(isfinite(xyzB_cache(:,1:2)), 2);
            finitePsi = isfinite(psi_cache);

            fprintf(['Projection status at iter %02d: totalOpt=%d | finiteAxy=%d | ' ...
                'finiteBxy=%d | finiteBothXY=%d | finitePsi=%d | lostPsi=%d\n'], ...
                a, numel(psi_cache), nnz(finiteA_xy), nnz(finiteB_xy), ...
                nnz(finiteA_xy & finiteB_xy), nnz(finitePsi), nnz(~finitePsi));

            metrics.debug.iterTotalOpt(a) = numel(psi_cache);
            metrics.debug.iterFiniteAxy(a) = nnz(finiteA_xy);
            metrics.debug.iterFiniteBxy(a) = nnz(finiteB_xy);
            metrics.debug.iterFiniteBothXY(a) = nnz(finiteA_xy & finiteB_xy);
            metrics.debug.iterFinitePsi(a) = nnz(finitePsi);
        end

        % ---------------------------------------------------------------------
        % Step 2: Build A-cell and B-cell lookups from current geometry
        % ---------------------------------------------------------------------
        base_idx_local = find(isfinite(psi_cache));
        if isempty(base_idx_local)
            metrics = KLT_update_convergence_metrics( ...
                metrics, ...
                meanAbsResidualHist, medianAbsResidualHist, rmsResidualHist, ...
                meanSignedResidualHist, medianSignedResidualHist, ...
                meanAbsResidualPostHist, medianAbsResidualPostHist, rmsResidualPostHist, ...
                meanSignedResidualPostHist, medianSignedResidualPostHist, ...
                medianAbsDzHist, maxAbsDzHist);

            checkpoint = KLT_build_checkpoint_Apoint( ...
                aa, idx, a, increment, alpha, maxIter, polarityMode, ...
                verbosePolarityPrint, jacobiMaxStep, jacobiMinSensitivity, ...
                jacobiRequireImprovement, jacobiFallbackToIncrement, ...
                phi, shift, e_s, e_n, origin, ...
                xi, yi, X_rot, Y_rot, Nx, Ny, dx, dy, x_order, y_order, ...
                active_window_mask, nSweepCells, opt_idx_full, ...
                cell_selected_idx_init, cell_median_val_init, wse_map, ...
                zi, zi_post, dzMapHist, signMapHist, residualMapHist, ...
                residualPostMapHist, selectedPathMapHist, xi_hist, yi_hist, ...
                meanAbsResidualHist, medianAbsResidualHist, rmsResidualHist, ...
                meanSignedResidualHist, medianSignedResidualHist, ...
                meanAbsResidualPostHist, medianAbsResidualPostHist, ...
                rmsResidualPostHist, meanSignedResidualPostHist, ...
                medianSignedResidualPostHist, medianAbsDzHist, maxAbsDzHist, ...
                activeCellFracHist, updatedCellFracHist, metrics);
            save(checkpointFile, 'checkpoint', '-v7.3');
            fprintf('Checkpoint saved: %s (completed iter %d)\n', checkpointFile, a);
            continue
        end

        abs_uv_for_local = opt_idx_full(base_idx_local);
        opt_uv_for_local = base_idx_local;

        xA_path = xyzA_cache(base_idx_local, 1);
        yA_path = xyzA_cache(base_idx_local, 2);
        xB_path = xyzB_cache(base_idx_local, 1);
        yB_path = xyzB_cache(base_idx_local, 2);
        r_path = psi_cache(base_idx_local);

        [cellA_x, cellA_y] = KLT_points_to_cells( ...
            xA_path, yA_path, xi, yi, dx, dy, origin, e_s, e_n);
        [cellB_x, cellB_y] = KLT_points_to_cells( ...
            xB_path, yB_path, xi, yi, dx, dy, origin, e_s, e_n);

        valid_bins = isfinite(cellA_x) & isfinite(cellA_y) & ...
            isfinite(cellB_x) & isfinite(cellB_y);

        finite_xy = isfinite(xA_path) & isfinite(yA_path) & ...
            isfinite(xB_path) & isfinite(yB_path);

        valid_geom = isfinite(r_path) & finite_xy & valid_bins;
        valid_diffcol = valid_geom & (cellA_x ~= cellB_x);

        lostToSameCol = valid_geom & (cellA_x == cellB_x);
        lostToOutOfGrid = isfinite(r_path) & finite_xy & ~valid_bins;

        if enforceCurrentDiffCol
            valid_path = valid_diffcol;
        else
            valid_path = valid_geom;
        end

        if debugDiagnostics
            fprintf(['Candidate status at iter %02d: geomValid=%d | diffColValid=%d | ' ...
                'lostToSameCol=%d | lostToOutOfGrid=%d | currentDiffCol=%d\n'], ...
                a, nnz(valid_geom), nnz(valid_diffcol), ...
                nnz(lostToSameCol), nnz(lostToOutOfGrid), enforceCurrentDiffCol);

            metrics.debug.iterGeomValid(a) = nnz(valid_geom);
            metrics.debug.iterDiffColValid(a) = nnz(valid_diffcol);
            metrics.debug.iterLostToSameCol(a) = nnz(lostToSameCol);
            metrics.debug.iterLostToOutOfGrid(a) = nnz(lostToOutOfGrid);
        end

        cellA_path_ids = cell(Ny, Nx);
        cellB_path_ids = cell(Ny, Nx);
        valid_idx = find(valid_path);

        if ~isempty(valid_idx)
            linA_ids = sub2ind([Ny, Nx], cellA_y(valid_idx), cellA_x(valid_idx));
            [linA_sorted, orderA] = sort(linA_ids);
            valid_idx_A_sorted = valid_idx(orderA);
            cutA = [1; find(diff(linA_sorted)) + 1; numel(linA_sorted) + 1];
            for g = 1:numel(cutA)-1
                i1 = cutA(g);
                i2 = cutA(g+1) - 1;
                cellA_path_ids{linA_sorted(i1)} = valid_idx_A_sorted(i1:i2).';
            end

            linB_ids = sub2ind([Ny, Nx], cellB_y(valid_idx), cellB_x(valid_idx));
            [linB_sorted, orderB] = sort(linB_ids);
            valid_idx_B_sorted = valid_idx(orderB);
            cutB = [1; find(diff(linB_sorted)) + 1; numel(linB_sorted) + 1];
            for g = 1:numel(cutB)-1
                i1 = cutB(g);
                i2 = cutB(g+1) - 1;
                cellB_path_ids{linB_sorted(i1)} = valid_idx_B_sorted(i1:i2).';
            end
        end

        if debugDiagnostics
            [candidateCountA_start, candidateCountB_start] = ...
                KLT_count_cell_path_ids(cellA_path_ids, cellB_path_ids);

            metrics.debug.candidateCountAStartHist{a} = candidateCountA_start;
            metrics.debug.candidateCountBStartHist{a} = candidateCountB_start;
        end

        active_rows_by_col = cell(1, Nx);
        nActive = 0;
        for ix_tmp = x_order
            mask_here = ~cellfun('isempty', cellA_path_ids(y_order, ix_tmp));
            rows_here = int32(y_order(mask_here));
            active_rows_by_col{ix_tmp} = rows_here(:);
            nActive = nActive + numel(rows_here);
        end

        fprintf('Active cells: %d / %d (%.1f%%)\n', ...
            nActive, nSweepCells, 100 * nActive / max(nSweepCells, 1));

        % ---------------------------------------------------------------------
        % Data-supported row window diagnostics
        % ---------------------------------------------------------------------
        rowActiveCount = zeros(Ny, 1);
        for ix_tmp = x_order
            rows_here_tmp = double(active_rows_by_col{ix_tmp});
            rowActiveCount(rows_here_tmp) = rowActiveCount(rows_here_tmp) + 1;
        end

        activeRowsAll = find(rowActiveCount > 0);

        if isempty(activeRowsAll)
            dataRowMin = NaN;
            dataRowMax = NaN;
            dataRows = [];
            nDataSupportedCells = NaN;
            lowBandRows = [];
            highBandRows = [];
            lowBandActive = 0;
            highBandActive = 0;
            fprintf('Data-supported row window: none\n');
        else
            dataRowMin = min(activeRowsAll);
            dataRowMax = max(activeRowsAll);
            dataRows = dataRowMin:dataRowMax;

            nDataSupportedCells = numel(dataRows) * numel(x_order);

            nLowUse = min(nNearDebugRows, numel(dataRows));
            lowBandRows = dataRows(1:nLowUse);
            highBandRows = dataRows(max(1, numel(dataRows)-nNearDebugRows+1):end);

            lowBandActive = sum(rowActiveCount(lowBandRows));
            highBandActive = sum(rowActiveCount(highBandRows));

            fprintf('Near-camera/data-start active cells in first %d active rows: %d\n', ...
                numel(lowBandRows), lowBandActive);

            fprintf('Active cells relative to full sweep:          %d / %d (%.1f%%)\n', ...
                nActive, nSweepCells, 100*nActive/max(nSweepCells,1));

            fprintf('Active cells relative to data-supported rows: %d / %d (%.1f%%)\n', ...
                nActive, nDataSupportedCells, 100*nActive/max(nDataSupportedCells,1));

            fprintf('Data-supported row window: rows %d to %d | rows=%d\n', ...
                dataRowMin, dataRowMax, numel(dataRows));

            fprintf('Active-low band rows:  %d to %d | activeCells=%d\n', ...
                lowBandRows(1), lowBandRows(end), lowBandActive);

            fprintf('Active-high band rows: %d to %d | activeCells=%d\n', ...
                highBandRows(1), highBandRows(end), highBandActive);
        end

        if debugDiagnostics
            metrics.debug.rowActiveCountStartHist{a} = ...
                uint16(min(rowActiveCount, double(intmax('uint16'))));

            if isempty(activeRowsAll)
                nearBandData = [];
                dataRowsContiguous = [];
                dataSupportedMask = false(Ny, Nx);
                nDataSupportedCellsDebug = 0;
                activeFracDataSupported = NaN;
            else
                nNearUseData = min(nNearDebugRows, numel(activeRowsAll));
                nearBandData = activeRowsAll(1:nNearUseData);
                dataRowsContiguous = dataRowMin:dataRowMax;
                dataSupportedMask = false(Ny, Nx);
                dataSupportedMask(dataRowsContiguous, x_order) = true;
                nDataSupportedCellsDebug = nnz(dataSupportedMask);
                activeFracDataSupported = nActive / max(nDataSupportedCellsDebug, 1);
            end

            nearActiveData = KLT_count_near_active_rows(active_rows_by_col, x_order, nearBandData);

            metrics.debug.nearBandDataHist{a} = nearBandData;
            metrics.debug.iterNearActiveDataCells(a) = nearActiveData;
            metrics.debug.dataRowMinHist(a) = dataRowMin;
            metrics.debug.dataRowMaxHist(a) = dataRowMax;
            metrics.debug.nDataSupportedCellsHist(a) = nDataSupportedCellsDebug;
            metrics.debug.activeFracDataSupportedHist(a) = activeFracDataSupported;
            metrics.debug.dataSupportedMaskHist{a} = dataSupportedMask;

            if ~isempty(activeRowsAll)
                nEdge = min(nNearDebugRows, numel(activeRowsAll));
                activeLowBand = activeRowsAll(1:nEdge);
                activeHighBand = activeRowsAll(end-nEdge+1:end);
                metrics.debug.activeLowBandHist{a} = activeLowBand;
                metrics.debug.activeHighBandHist{a} = activeHighBand;
            end

            nearActive = KLT_count_near_active_rows(active_rows_by_col, x_order, nearBand);
            metrics.debug.iterActiveCells(a) = nActive;
            metrics.debug.iterNearActiveCells(a) = nearActive;

            fprintf('Near-camera active cells in first %d sweep rows: %d\n', ...
                numel(nearBand), nearActive);
        end

        % ---------------------------------------------------------------------
        % Iteration diagnostics
        % ---------------------------------------------------------------------
        zi{a} = nan(Ny, Nx);
        zi_post{a} = nan(Ny, Nx);
        dz_iter = nan(Ny, Nx);
        sign_iter = nan(Ny, Nx);
        selected_path_iter = nan(Ny, Nx);
        updated_mask_iter = false(Ny, Nx);
        eligible_mask_iter = false(Ny, Nx);
        eligibleNoUpdate_iter = false(Ny, Nx);
        reject_code_iter = nan(Ny, Nx);

        residual_used_pre = nan(nSweepCells, 1);
        residual_used_post = nan(nSweepCells, 1);
        nResidualUsed = 0;

        selectedOptIdxThisIter = [];

        % =====================================================================
        % Column-wise Jacobi updates, swept downstream column by column
        % =====================================================================
        for ix_now = x_order
            rows_now = active_rows_by_col{ix_now};
            if isempty(rows_now)
                continue
            end

            dem_col_frozen = wse_map{aa,idx};
            r_path_frozen = r_path;
            nRows = numel(rows_now);

            % Preselect one representative A-assigned path per active row
            pre_eligible = false(nRows, 1);
            pre_p_sel = nan(nRows, 1);
            pre_opt_uv_idx = nan(nRows, 1);
            pre_abs_uv_idx = nan(nRows, 1);
            pre_r0 = nan(nRows, 1);
            pre_uvA = nan(nRows, size(xyzA_opt_uv, 2));
            pre_uvB = nan(nRows, size(xyzB_opt_uv, 2));

            for j = 1:nRows
                [pre_eligible(j), ...
                    pre_p_sel(j), ...
                    pre_opt_uv_idx(j), ...
                    pre_abs_uv_idx(j), ...
                    pre_r0(j), ...
                    pre_uvA(j,:), ...
                    pre_uvB(j,:)] = KLT_preselect_row_candidate_scalar( ...
                    cellA_path_ids{rows_now(j), ix_now}, ...
                    r_path_frozen, abs_uv_for_local, opt_uv_for_local, ...
                    xyzA_opt_uv, xyzB_opt_uv);
            end

            camA_local = app.camA;
            X_rot_local = X_rot;
            Y_rot_local = Y_rot;
            phi_local = phi;
            shift_local = shift;
            flowAngleOffsetDeg_local = flowAngleOffsetDeg;
            increment_local = increment;
            usePreparedDEMRowTrials_local = usePreparedDEMRowTrials;
            demPrepCurrent_local = demPrepCurrent;
            jacobiMaxStep_local = jacobiMaxStep;
            jacobiMinSensitivity_local = jacobiMinSensitivity;
            jacobiRequireImprovement_local = jacobiRequireImprovement;
            jacobiFallbackToIncrement_local = jacobiFallbackToIncrement;

            eligible_col = false(nRows, 1);
            psel_col = nan(nRows, 1);
            opt_uv_idx_col = nan(nRows, 1);
            abs_uv_idx_col = nan(nRows, 1);
            r0_col = nan(nRows, 1);
            dz_req_col = zeros(nRows, 1);
            r_plus_col = nan(nRows, 1);
            r_minus_col = nan(nRows, 1);
            r_solve_col = nan(nRows, 1);
            sens_col = nan(nRows, 1);
            trialCode_col = zeros(nRows, 1, 'uint8');
            rejectCode_col = zeros(nRows, 1, 'uint8');

            if useParallelRows
                parfor j = 1:nRows
                    [eligible_col(j), ...
                        psel_col(j), ...
                        opt_uv_idx_col(j), ...
                        abs_uv_idx_col(j), ...
                        r0_col(j), ...
                        dz_req_col(j), ...
                        r_plus_col(j), ...
                        r_minus_col(j), ...
                        r_solve_col(j), ...
                        sens_col(j), ...
                        trialCode_col(j), ...
                        rejectCode_col(j)] = KLT_evaluate_preselected_row_candidate_jacobi_scalar( ...
                        pre_eligible(j), ...
                        pre_p_sel(j), ...
                        pre_opt_uv_idx(j), ...
                        pre_abs_uv_idx(j), ...
                        pre_r0(j), ...
                        pre_uvA(j,:), ...
                        pre_uvB(j,:), ...
                        camA_local, dem_col_frozen, ...
                        usePreparedDEMRowTrials_local, demPrepCurrent_local, ...
                        ix_now, rows_now(j), increment_local, ...
                        jacobiMaxStep_local, jacobiMinSensitivity_local, ...
                        jacobiRequireImprovement_local, jacobiFallbackToIncrement_local, ...
                        X_rot_local, Y_rot_local, ...
                        phi_local, shift_local, flowAngleOffsetDeg_local);
                end
            else
                for j = 1:nRows
                    [eligible_col(j), ...
                        psel_col(j), ...
                        opt_uv_idx_col(j), ...
                        abs_uv_idx_col(j), ...
                        r0_col(j), ...
                        dz_req_col(j), ...
                        r_plus_col(j), ...
                        r_minus_col(j), ...
                        r_solve_col(j), ...
                        sens_col(j), ...
                        trialCode_col(j), ...
                        rejectCode_col(j)] = KLT_evaluate_preselected_row_candidate_jacobi_scalar( ...
                        pre_eligible(j), ...
                        pre_p_sel(j), ...
                        pre_opt_uv_idx(j), ...
                        pre_abs_uv_idx(j), ...
                        pre_r0(j), ...
                        pre_uvA(j,:), ...
                        pre_uvB(j,:), ...
                        camA_local, dem_col_frozen, ...
                        usePreparedDEMRowTrials_local, demPrepCurrent_local, ...
                        ix_now, rows_now(j), increment_local, ...
                        jacobiMaxStep_local, jacobiMinSensitivity_local, ...
                        jacobiRequireImprovement_local, jacobiFallbackToIncrement_local, ...
                        X_rot_local, Y_rot_local, ...
                        phi_local, shift_local, flowAngleOffsetDeg_local);
                end
            end

            updated_col = eligible_col & isfinite(dz_req_col) & dz_req_col ~= 0;

            goodRep = eligible_col & isfinite(opt_uv_idx_col);
            if any(goodRep)
                selectedOptIdxThisIter = [selectedOptIdxThisIter; opt_uv_idx_col(goodRep)]; %#ok<AGROW>
            end

            if verbosePolarityPrint
                for j = 1:nRows
                    if ~eligible_col(j) || ~isfinite(psel_col(j))
                        continue
                    end

                    iy_now = rows_now(j);
                    abs_uv_idx = abs_uv_idx_col(j);
                    r0 = r0_col(j);
                    r_plus = r_plus_col(j);
                    r_minus = r_minus_col(j);
                    r_solve = r_solve_col(j);
                    sens_here = sens_col(j);
                    dz_req = dz_req_col(j);

                    if dz_req > 0
                        chosenLabel = '+dz';
                    elseif dz_req < 0
                        chosenLabel = '-dz';
                    else
                        chosenLabel = '0';
                    end

                    trialOrderText = KLT_jacobi_trial_code_to_text(trialCode_col(j));
                    rejectText = KLT_reject_code_to_text(rejectCode_col(j));

                    fprintf( ...
                        ['Iter %02d | cell(row=%3d,col=%3d) | uvIdx=%6d | ' ...
                        'r0=%+10.6f | r(+)=%+10.6f | r(-)=%+10.6f | ' ...
                        'drdz=%+12.6e | dzReq=%+10.6f | rSolve=%+10.6f | ' ...
                        'mode=%s | rejectCode=%d:%s | choose %s\n'], ...
                        a, iy_now, ix_now, abs_uv_idx, ...
                        r0, r_plus, r_minus, sens_here, dz_req, r_solve, ...
                        trialOrderText, rejectCode_col(j), rejectText, chosenLabel);
                end
            end

            % Apply all accepted row updates for this column together
            rows_apply = rows_now(updated_col);
            dz_apply = alpha * dz_req_col(updated_col);

            if ~isempty(rows_apply)
                lin_apply = sub2ind([Ny, Nx], double(rows_apply), ...
                    repmat(ix_now, numel(rows_apply), 1));
                wse_map{aa,idx}(lin_apply) = ...
                    wse_map{aa,idx}(lin_apply) + dz_apply;

                if usePreparedDEMSeed || usePreparedDEMRefresh
                    demPrepCurrent = app.camA.updatePreparedDEM(demPrepBase, wse_map{aa,idx});
                end
            end

            for j = 1:nRows
                iy_now = rows_now(j);
                if eligible_col(j)
                    eligible_mask_iter(iy_now, ix_now) = true;
                    reject_code_iter(iy_now, ix_now) = double(rejectCode_col(j));
                end
            end

            % Affected-cache refresh using both A-cell and B-cell membership
            if ~isempty(rows_apply)
                if useStencilAffectedRefresh
                    stencil_lin = [];
                    for q = 1:numel(rows_apply)
                        rr0 = double(rows_apply(q));
                        cc0 = ix_now;
                        rr = max(1, rr0-1):min(Ny, rr0+1);
                        cc = max(1, cc0-1):min(Nx, cc0+1);
                        [CC, RR] = meshgrid(cc, rr);
                        stencil_lin = [stencil_lin; sub2ind([Ny, Nx], RR(:), CC(:))]; %#ok<AGROW>
                    end
                    affected_cell_lin = unique(stencil_lin, 'stable');
                else
                    affected_cell_lin = sub2ind([Ny, Nx], double(rows_apply(:)), ...
                        repmat(ix_now, numel(rows_apply), 1));
                end

                tmp_cells = cell(2*numel(affected_cell_lin), 1);
                for q = 1:numel(affected_cell_lin)
                    tmp_cells{2*q - 1} = cellA_path_ids{affected_cell_lin(q)}(:);
                    tmp_cells{2*q} = cellB_path_ids{affected_cell_lin(q)}(:);
                end

                affected_local_col = vertcat(tmp_cells{:});
                if ~isempty(affected_local_col)
                    affected_local_col = unique(affected_local_col, 'stable');
                    affected_cache_col = base_idx_local(affected_local_col);

                    if usePreparedDEMRefresh
                        if useWarmStartsRefresh
                            xy0A_aff = xyzA_cache(affected_cache_col, 1:2);
                            xy0B_aff = xyzB_cache(affected_cache_col, 1:2);
                            [xyzA_aff, xyzB_aff] = KLT_invproject_uv_pairs_safe_or_prepared( ...
                                app.camA, ...
                                xyzA_opt_uv(affected_cache_col, :), ...
                                xyzB_opt_uv(affected_cache_col, :), ...
                                X_rot, Y_rot, wse_map{aa,idx}, ...
                                true, demPrepCurrent, ...
                                xy0A_aff, xy0B_aff);
                        else
                            [xyzA_aff, xyzB_aff] = KLT_invproject_uv_pairs_safe_or_prepared( ...
                                app.camA, ...
                                xyzA_opt_uv(affected_cache_col, :), ...
                                xyzB_opt_uv(affected_cache_col, :), ...
                                X_rot, Y_rot, wse_map{aa,idx}, ...
                                true, demPrepCurrent);
                        end
                    else
                        [xyzA_aff, xyzB_aff] = KLT_invproject_uv_pairs_safe_or_prepared( ...
                            app.camA, ...
                            xyzA_opt_uv(affected_cache_col, :), ...
                            xyzB_opt_uv(affected_cache_col, :), ...
                            X_rot, Y_rot, wse_map{aa,idx}, ...
                            false, demPrepCurrent);
                    end

                    xyzA_cache(affected_cache_col, :) = xyzA_aff;
                    xyzB_cache(affected_cache_col, :) = xyzB_aff;

                    psi_cache(affected_cache_col) = KLT_compute_residual_from_projected( ...
                        xyzA_cache(affected_cache_col, 1:2), ...
                        xyzB_cache(affected_cache_col, 1:2), phi, shift, flowAngleOffsetDeg);
                    r_path(affected_local_col) = psi_cache(affected_cache_col);
                end
            end

            % Record diagnostics for all eligible cells
            for j = 1:nRows
                if ~eligible_col(j) || ~isfinite(psel_col(j))
                    continue
                end

                iy_now = rows_now(j);
                p_sel = psel_col(j);
                abs_uv_idx = abs_uv_idx_col(j);
                r0 = r0_col(j);
                r_true_post = r_path(p_sel);

                zi{a}(iy_now, ix_now) = r0;
                zi_post{a}(iy_now, ix_now) = r_true_post;
                selected_path_iter(iy_now, ix_now) = abs_uv_idx;

                if updated_col(j)
                    dz_iter(iy_now, ix_now) = alpha * dz_req_col(j);
                    sign_iter(iy_now, ix_now) = sign(dz_req_col(j));
                    updated_mask_iter(iy_now, ix_now) = true;
                else
                    eligibleNoUpdate_iter(iy_now, ix_now) = true;
                end

                nResidualUsed = nResidualUsed + 1;
                residual_used_pre(nResidualUsed, 1) = r0;
                residual_used_post(nResidualUsed, 1) = r_true_post;
            end

            nEligibleCol = nnz(eligible_col);
            nUpdatedCol = nnz(updated_col);
            rows_eligible = double(rows_now(eligible_col));

            if debugDiagnostics
                near_col_mask = ismember(double(rows_now), nearBand(:));
                nearEligibleCol = nnz(eligible_col & near_col_mask);
                nearUpdatedCol = nnz(updated_col & near_col_mask);
            else
                nearEligibleCol = NaN;
                nearUpdatedCol = NaN;
            end

            if ~isempty(rows_eligible)
                meanAbsPreCol = mean(abs(r0_col(eligible_col)), 'omitnan');
                lin_eligible = sub2ind([Ny, Nx], rows_eligible(:), ...
                    repmat(ix_now, numel(rows_eligible), 1));
                meanAbsPostCol = mean(abs(zi_post{a}(lin_eligible)), 'omitnan');
            else
                meanAbsPreCol = NaN;
                meanAbsPostCol = NaN;
            end

            fprintf(['Iter %02d | col=%4d | eligibleRows=%4d | updatedRows=%4d | ' ...
                'mean|r| pre=%10.6f | mean|r| post=%10.6f\n'], ...
                a, ix_now, nEligibleCol, nUpdatedCol, meanAbsPreCol, meanAbsPostCol);

            if debugDiagnostics && nearEligibleCol > 0
                fprintf('          nearRows: eligible=%4d | updated=%4d\n', ...
                    nearEligibleCol, nearUpdatedCol);
            end

            % Dynamic membership rebuild for future columns
            if rebuildMembershipAfterChangedColumn && ~isempty(rows_apply)
                [base_idx_local, ...
                    abs_uv_for_local, ...
                    opt_uv_for_local, ...
                    r_path, ...
                    cellA_path_ids, ...
                    cellB_path_ids, ...
                    active_rows_by_col, ...
                    nActive] = KLT_rebuild_AB_membership_from_cache( ...
                    xyzA_cache, xyzB_cache, psi_cache, opt_idx_full, ...
                    xi, yi, dx, dy, origin, e_s, e_n, ...
                    Ny, Nx, x_order, y_order, enforceCurrentDiffCol);

                fprintf('Membership rebuilt after col=%4d | activeCells now=%d / %d (%.1f%%)\n', ...
                    ix_now, nActive, nSweepCells, 100 * nActive / max(nSweepCells, 1));
            end
        end

        % ---------------------------------------------------------------------
        % Store iteration-level maps and standard residual/update metrics
        % ---------------------------------------------------------------------
        residual_used_pre = residual_used_pre(1:nResidualUsed);
        residual_used_post = residual_used_post(1:nResidualUsed);

        residualMapHist{a} = zi{a};
        residualPostMapHist{a} = zi_post{a};
        dzMapHist{a} = dz_iter;
        signMapHist{a} = sign_iter;
        selectedPathMapHist{a} = selected_path_iter;
        xi_hist{a} = xi;
        yi_hist{a} = yi;

        if debugDiagnostics
            [candidateCountA_end, candidateCountB_end] = ...
                KLT_count_cell_path_ids(cellA_path_ids, cellB_path_ids);

            metrics.debug.candidateCountAEndHist{a} = candidateCountA_end;
            metrics.debug.candidateCountBEndHist{a} = candidateCountB_end;
            metrics.debug.eligibleNoUpdateHist{a} = eligibleNoUpdate_iter;
        end

        if ~isempty(residual_used_pre)
            meanAbsResidualHist(a) = mean(abs(residual_used_pre), 'omitnan');
            medianAbsResidualHist(a) = median(abs(residual_used_pre), 'omitnan');
            rmsResidualHist(a) = sqrt(mean(residual_used_pre.^2, 'omitnan'));
            meanSignedResidualHist(a) = mean(residual_used_pre, 'omitnan');
            medianSignedResidualHist(a) = median(residual_used_pre, 'omitnan');
        end

        if ~isempty(residual_used_post)
            meanAbsResidualPostHist(a) = mean(abs(residual_used_post), 'omitnan');
            medianAbsResidualPostHist(a) = median(abs(residual_used_post), 'omitnan');
            rmsResidualPostHist(a) = sqrt(mean(residual_used_post.^2, 'omitnan'));
            meanSignedResidualPostHist(a) = mean(residual_used_post, 'omitnan');
            medianSignedResidualPostHist(a) = median(residual_used_post, 'omitnan');
        end

        dzVals = dz_iter(isfinite(dz_iter) & dz_iter ~= 0);
        if ~isempty(dzVals)
            medianAbsDzHist(a) = median(abs(dzVals), 'omitnan');
            maxAbsDzHist(a) = max(abs(dzVals));
        end

        if nSweepCells > 0
            activeCellFracHist(a) = ...
                nnz(eligible_mask_iter & active_window_mask) / nSweepCells;
            updatedCellFracHist(a) = ...
                nnz(updated_mask_iter & active_window_mask) / nSweepCells;
        end

        if debugDiagnostics
            metrics.debug.iterEligibleCells(a) = nnz(eligible_mask_iter & active_window_mask);
            metrics.debug.iterUpdatedCells(a) = nnz(updated_mask_iter & active_window_mask);
        end

        % Compact iteration summary
        if ~isempty(dataRows)
            dataWindowMask = false(Ny, Nx);
            dataWindowMask(dataRows, x_order) = true;

            lowBandMask = false(Ny, Nx);
            lowBandMask(lowBandRows, x_order) = true;

            highBandMask = false(Ny, Nx);
            highBandMask(highBandRows, x_order) = true;

            nEligibleData = nnz(eligible_mask_iter & dataWindowMask);
            nUpdatedData = nnz(updated_mask_iter & dataWindowMask);

            nEligibleLow = nnz(eligible_mask_iter & lowBandMask);
            nUpdatedLow = nnz(updated_mask_iter & lowBandMask);

            nEligibleHigh = nnz(eligible_mask_iter & highBandMask);
            nUpdatedHigh = nnz(updated_mask_iter & highBandMask);

            fprintf(['Iter %02d summary | dataRows=%d:%d | dataActive=%d/%d %.1f%% | ' ...
                'lowBand active=%d eligible=%d updated=%d | ' ...
                'highBand active=%d eligible=%d updated=%d | ' ...
                'allData eligible=%d updated=%d\n'], ...
                a, dataRowMin, dataRowMax, ...
                nActive, nDataSupportedCells, 100*nActive/max(nDataSupportedCells,1), ...
                lowBandActive, nEligibleLow, nUpdatedLow, ...
                highBandActive, nEligibleHigh, nUpdatedHigh, ...
                nEligibleData, nUpdatedData);

            if debugDiagnostics
                metrics.debug.iterDataRowMin(a) = dataRowMin;
                metrics.debug.iterDataRowMax(a) = dataRowMax;
                metrics.debug.iterDataSupportedCells(a) = nDataSupportedCells;
                metrics.debug.iterLowBandActive(a) = lowBandActive;
                metrics.debug.iterHighBandActive(a) = highBandActive;
                metrics.debug.iterLowBandEligible(a) = nEligibleLow;
                metrics.debug.iterLowBandUpdated(a) = nUpdatedLow;
                metrics.debug.iterHighBandEligible(a) = nEligibleHigh;
                metrics.debug.iterHighBandUpdated(a) = nUpdatedHigh;
            end
        end

        % ---------------------------------------------------------------------
        % New advanced convergence diagnostics
        % ---------------------------------------------------------------------
        if doEndIterationFullReprojection
            [xyzA_end, xyzB_end] = KLT_invproject_uv_pairs_safe_or_prepared( ...
                app.camA, ...
                xyzA_opt_uv, xyzB_opt_uv, ...
                X_rot, Y_rot, wse_map{aa,idx}, ...
                usePreparedDEMSeed, demPrepCurrent);

            psi_end_all = KLT_compute_residual_from_projected( ...
                xyzA_end(:,1:2), xyzB_end(:,1:2), ...
                phi, shift, flowAngleOffsetDeg);
        else
            xyzA_end = xyzA_cache;
            xyzB_end = xyzB_cache;
            psi_end_all = psi_cache;
        end

        [endRepResidualMap, endRepSelectedLocalIdx, endValidPathMask] = ...
            KLT_build_end_iteration_representative_map( ...
            xyzA_end, xyzB_end, psi_end_all, ...
            xi, yi, dx, dy, origin, e_s, e_n, ...
            Ny, Nx, enforceCurrentDiffCol);

        metrics = KLT_update_advanced_convergence_diagnostics( ...
            metrics, a, maxIter, ...
            residual_used_pre, residual_used_post, ...
            dz_iter, dzMapHist, wse_map{aa,idx-1}, wse_map{aa,idx}, ...
            nSweepCells, eligible_mask_iter & active_window_mask, ...
            y_order, x_order, Ny, Nx, ...
            psi_end_all, endRepResidualMap, endRepSelectedLocalIdx, endValidPathMask, ...
            selectedOptIdxThisIter, reject_code_iter, ...
            nConvergenceBands, convergenceWindow);

        % Warm starts from full end-of-iteration projection
        prev_xyA_opt = xyzA_end(:,1:2);
        prev_xyB_opt = xyzB_end(:,1:2);

        metrics = KLT_update_convergence_metrics( ...
            metrics, ...
            meanAbsResidualHist, medianAbsResidualHist, rmsResidualHist, ...
            meanSignedResidualHist, medianSignedResidualHist, ...
            meanAbsResidualPostHist, medianAbsResidualPostHist, rmsResidualPostHist, ...
            meanSignedResidualPostHist, medianSignedResidualPostHist, ...
            medianAbsDzHist, maxAbsDzHist);

        checkpoint = KLT_build_checkpoint_Apoint( ...
            aa, idx, a, increment, alpha, maxIter, polarityMode, ...
            verbosePolarityPrint, jacobiMaxStep, jacobiMinSensitivity, ...
            jacobiRequireImprovement, jacobiFallbackToIncrement, ...
            phi, shift, e_s, e_n, origin, ...
            xi, yi, X_rot, Y_rot, Nx, Ny, dx, dy, x_order, y_order, ...
            active_window_mask, nSweepCells, opt_idx_full, ...
            cell_selected_idx_init, cell_median_val_init, wse_map, ...
            zi, zi_post, dzMapHist, signMapHist, residualMapHist, ...
            residualPostMapHist, selectedPathMapHist, xi_hist, yi_hist, ...
            meanAbsResidualHist, medianAbsResidualHist, rmsResidualHist, ...
            meanSignedResidualHist, medianSignedResidualHist, ...
            meanAbsResidualPostHist, medianAbsResidualPostHist, ...
            rmsResidualPostHist, meanSignedResidualPostHist, ...
            medianSignedResidualPostHist, medianAbsDzHist, maxAbsDzHist, ...
            activeCellFracHist, updatedCellFracHist, metrics);
        save(checkpointFile, 'checkpoint', '-v7.3');
        fprintf('Checkpoint saved: %s (completed iter %d)\n', checkpointFile, a);
    end

    % -------------------------------------------------------------------------
    % Final full reprojection on final surface
    % -------------------------------------------------------------------------
    metrics.xyzA_all_final = app.camA.invproject( ...
        xyzA_wse{aa}, X_rot, Y_rot, wse_map{aa,idx});
    metrics.xyzB_all_final = app.camA.invproject( ...
        xyzB_wse{aa}, X_rot, Y_rot, wse_map{aa,idx});

    metrics.residualMapLast = zi{KLT_last_valid_iter(zi)};
    metrics.residualPostMapLast = zi_post{KLT_last_valid_iter(zi_post)};
    metrics.dzMapLast = dzMapHist{KLT_last_valid_iter(dzMapHist)};
    metrics.signMapLast = signMapHist{KLT_last_valid_iter(signMapHist)};
    metrics.selectedPathMapLast = selectedPathMapHist{KLT_last_valid_iter(selectedPathMapHist)};

    % Final metrics bundle
    metrics.increment = increment;
    metrics.alpha = alpha;
    metrics.maxIter = maxIter;
    metrics.jacobiMaxStep = jacobiMaxStep;
    metrics.jacobiMinSensitivity = jacobiMinSensitivity;
    metrics.jacobiRequireImprovement = jacobiRequireImprovement;
    metrics.jacobiFallbackToIncrement = jacobiFallbackToIncrement;
    metrics.flowAngleOffsetDeg = flowAngleOffsetDeg;
    metrics.flowAngleOffsetDegWasInput = flowAngleOffsetDegWasInput;
    metrics.flowAngleOffsetWindowSource = flowAngleOffsetWindowSource;
    metrics.flowAngleOffsetMode = flowAngleOffsetMode;
    metrics.residualZeroRowRange = residualZeroRowRange;
    metrics.residualZeroColRange = residualZeroColRange;
    metrics.residualZeroMedianRawDeg = residualZeroMedianRawDeg;
    metrics.activeCellFracHist = activeCellFracHist;
    metrics.updatedCellFracHist = updatedCellFracHist;
    metrics.residualMapHist = residualMapHist;
    metrics.residualPostMapHist = residualPostMapHist;
    metrics.dzMapHist = dzMapHist;
    metrics.signMapHist = signMapHist;
    metrics.selectedPathMapHist = selectedPathMapHist;
    metrics.xiHist = xi_hist;
    metrics.yiHist = yi_hist;

    metrics = KLT_update_convergence_metrics( ...
        metrics, ...
        meanAbsResidualHist, medianAbsResidualHist, rmsResidualHist, ...
        meanSignedResidualHist, medianSignedResidualHist, ...
        meanAbsResidualPostHist, medianAbsResidualPostHist, rmsResidualPostHist, ...
        meanSignedResidualPostHist, medianSignedResidualPostHist, ...
        medianAbsDzHist, maxAbsDzHist);

    metrics.activeCellFracFinal = KLT_last_finite(activeCellFracHist);
    metrics.updatedCellFracFinal = KLT_last_finite(updatedCellFracHist);
    metrics.activeCellFracAllIters = mean(activeCellFracHist, 'omitnan');
    metrics.updatedCellFracAllIters = mean(updatedCellFracHist, 'omitnan');

    checkpoint = KLT_build_checkpoint_Apoint( ...
        aa, idx, maxIter, increment, alpha, maxIter, polarityMode, ...
        verbosePolarityPrint, jacobiMaxStep, jacobiMinSensitivity, ...
        jacobiRequireImprovement, jacobiFallbackToIncrement, ...
        phi, shift, e_s, e_n, origin, ...
        xi, yi, X_rot, Y_rot, Nx, Ny, dx, dy, x_order, y_order, ...
        active_window_mask, nSweepCells, opt_idx_full, ...
        cell_selected_idx_init, cell_median_val_init, wse_map, ...
        zi, zi_post, dzMapHist, signMapHist, residualMapHist, ...
        residualPostMapHist, selectedPathMapHist, xi_hist, yi_hist, ...
        meanAbsResidualHist, medianAbsResidualHist, rmsResidualHist, ...
        meanSignedResidualHist, medianSignedResidualHist, ...
        meanAbsResidualPostHist, medianAbsResidualPostHist, ...
        rmsResidualPostHist, meanSignedResidualPostHist, ...
        medianSignedResidualPostHist, medianAbsDzHist, maxAbsDzHist, ...
        activeCellFracHist, updatedCellFracHist, metrics);
    checkpoint.nextIter = maxIter + 1;
    save(checkpointFile, 'checkpoint', '-v7.3');
    fprintf('Final checkpoint saved: %s\n', checkpointFile);
end

% =========================================================================
% Detect whether parfor can be used safely on this MATLAB install
% =========================================================================
function tf = KLT_can_use_parallel()
    tf = ~isempty(ver('parallel')) && license('test', 'Distrib_Computing_Toolbox');
end

% =========================================================================
% Checkpoint convenience helpers
% =========================================================================
function C = KLT_checkpoint_cell_or_default(checkpoint, fieldName, maxIter)
    if isfield(checkpoint, fieldName) && iscell(checkpoint.(fieldName))
        C = checkpoint.(fieldName);
        C = C(:);
        if numel(C) < maxIter
            C(end+1:maxIter,1) = {[]};
        elseif numel(C) > maxIter
            C = C(1:maxIter);
        end
    else
        C = cell(maxIter,1);
    end
end

function v = KLT_checkpoint_vec_or_default(checkpoint, fieldName, maxIter)
    v = nan(maxIter,1);
    if isfield(checkpoint, fieldName) && isnumeric(checkpoint.(fieldName))
        old = checkpoint.(fieldName);
        old = old(:);
        nCopy = min(numel(old), maxIter);
        v(1:nCopy) = old(1:nCopy);
    end
end

% =========================================================================
% Build checkpoint struct
% =========================================================================
function checkpoint = KLT_build_checkpoint_Apoint( ...
    aa, idx, iterCompleted, increment, alpha, maxIter, polarityMode, ...
    verbosePolarityPrint, jacobiMaxStep, jacobiMinSensitivity, ...
    jacobiRequireImprovement, jacobiFallbackToIncrement, ...
    phi, shift, e_s, e_n, origin, ...
    xi, yi, X_rot, Y_rot, Nx, Ny, dx, dy, x_order, y_order, ...
    active_window_mask, nSweepCells, opt_idx_full, ...
    cell_selected_idx_init, cell_median_val_init, wse_map, ...
    zi, zi_post, dzMapHist, signMapHist, residualMapHist, ...
    residualPostMapHist, selectedPathMapHist, xi_hist, yi_hist, ...
    meanAbsResidualHist, medianAbsResidualHist, rmsResidualHist, ...
    meanSignedResidualHist, medianSignedResidualHist, ...
    meanAbsResidualPostHist, medianAbsResidualPostHist, ...
    rmsResidualPostHist, meanSignedResidualPostHist, ...
    medianSignedResidualPostHist, medianAbsDzHist, maxAbsDzHist, ...
    activeCellFracHist, updatedCellFracHist, metrics_static)

    checkpoint = struct;
    checkpoint.aa = aa;
    checkpoint.idx = idx;
    checkpoint.lastCompletedIter = iterCompleted;
    checkpoint.nextIter = iterCompleted + 1;
    checkpoint.increment = increment;
    checkpoint.alpha = alpha;
    checkpoint.maxIter = maxIter;
    checkpoint.polarityMode = polarityMode;
    checkpoint.verbosePolarityPrint = verbosePolarityPrint;
    checkpoint.jacobiMaxStep = jacobiMaxStep;
    checkpoint.jacobiMinSensitivity = jacobiMinSensitivity;
    checkpoint.jacobiRequireImprovement = jacobiRequireImprovement;
    checkpoint.jacobiFallbackToIncrement = jacobiFallbackToIncrement;
    checkpoint.phi = phi;
    checkpoint.shift = shift;
    checkpoint.e_s = e_s;
    checkpoint.e_n = e_n;
    checkpoint.origin = origin;
    checkpoint.xi = xi;
    checkpoint.yi = yi;
    checkpoint.X_rot = X_rot;
    checkpoint.Y_rot = Y_rot;
    checkpoint.Nx = Nx;
    checkpoint.Ny = Ny;
    checkpoint.dx = dx;
    checkpoint.dy = dy;
    checkpoint.x_order = x_order;
    checkpoint.y_order = y_order;
    checkpoint.active_window_mask = active_window_mask;
    checkpoint.nSweepCells = nSweepCells;
    checkpoint.opt_idx_full = opt_idx_full;
    checkpoint.cell_selected_idx_init = cell_selected_idx_init;
    checkpoint.cell_median_val_init = cell_median_val_init;
    checkpoint.wse_map = wse_map;
    checkpoint.zi = zi;
    checkpoint.zi_post = zi_post;
    checkpoint.dzMapHist = dzMapHist;
    checkpoint.signMapHist = signMapHist;
    checkpoint.residualMapHist = residualMapHist;
    checkpoint.residualPostMapHist = residualPostMapHist;
    checkpoint.selectedPathMapHist = selectedPathMapHist;
    checkpoint.xi_hist = xi_hist;
    checkpoint.yi_hist = yi_hist;
    checkpoint.meanAbsResidualHist = meanAbsResidualHist;
    checkpoint.medianAbsResidualHist = medianAbsResidualHist;
    checkpoint.rmsResidualHist = rmsResidualHist;
    checkpoint.meanSignedResidualHist = meanSignedResidualHist;
    checkpoint.medianSignedResidualHist = medianSignedResidualHist;
    checkpoint.meanAbsResidualPostHist = meanAbsResidualPostHist;
    checkpoint.medianAbsResidualPostHist = medianAbsResidualPostHist;
    checkpoint.rmsResidualPostHist = rmsResidualPostHist;
    checkpoint.meanSignedResidualPostHist = meanSignedResidualPostHist;
    checkpoint.medianSignedResidualPostHist = medianSignedResidualPostHist;
    checkpoint.medianAbsDzHist = medianAbsDzHist;
    checkpoint.maxAbsDzHist = maxAbsDzHist;
    checkpoint.activeCellFracHist = activeCellFracHist;
    checkpoint.updatedCellFracHist = updatedCellFracHist;
    checkpoint.metrics_static = metrics_static;

    if isfield(metrics_static, 'flowAngleOffsetDeg')
        checkpoint.flowAngleOffsetDeg = metrics_static.flowAngleOffsetDeg;
    else
        checkpoint.flowAngleOffsetDeg = 0;
    end

    if isfield(metrics_static, 'residualZeroRowRange')
        checkpoint.residualZeroRowRange = metrics_static.residualZeroRowRange;
    else
        checkpoint.residualZeroRowRange = [];
    end

    if isfield(metrics_static, 'residualZeroColRange')
        checkpoint.residualZeroColRange = metrics_static.residualZeroColRange;
    else
        checkpoint.residualZeroColRange = [];
    end

    if isfield(metrics_static, 'flowAngleOffsetWindowSource')
        checkpoint.flowAngleOffsetWindowSource = metrics_static.flowAngleOffsetWindowSource;
    else
        checkpoint.flowAngleOffsetWindowSource = 'full_domain';
    end

    if isfield(metrics_static, 'flowAngleOffsetDegWasInput')
        checkpoint.flowAngleOffsetDegWasInput = metrics_static.flowAngleOffsetDegWasInput;
    else
        checkpoint.flowAngleOffsetDegWasInput = false;
    end

    if isfield(metrics_static, 'flowAngleOffsetMode')
        checkpoint.flowAngleOffsetMode = metrics_static.flowAngleOffsetMode;
    else
        checkpoint.flowAngleOffsetMode = 'residual_median_inferred';
    end

    if isfield(metrics_static, 'residualZeroMedianRawDeg')
        checkpoint.residualZeroMedianRawDeg = metrics_static.residualZeroMedianRawDeg;
    else
        checkpoint.residualZeroMedianRawDeg = NaN;
    end
end

% =========================================================================
% Evaluate one preselected row candidate using a column-wise Jacobi/Newton step
% =========================================================================
function [eligible, p_sel, opt_uv_idx, abs_uv_idx, r0, dz_required, ...
    r_plus, r_minus, r_solve, sensitivity, trialCode, rejectCode] = ...
    KLT_evaluate_preselected_row_candidate_jacobi_scalar( ...
    pre_eligible, pre_p_sel, pre_opt_uv_idx, pre_abs_uv_idx, pre_r0, ...
    pre_uvA, pre_uvB, ...
    camA, dem_current, usePreparedDEMRowTrials, demPrepCurrent, ...
    ix_now, iy_now, increment, ...
    jacobiMaxStep, jacobiMinSensitivity, ...
    jacobiRequireImprovement, jacobiFallbackToIncrement, ...
    TransX, TransY, phi, shift, flowAngleOffsetDeg)

    eligible = false;
    p_sel = NaN;
    opt_uv_idx = NaN;
    abs_uv_idx = NaN;
    r0 = NaN;
    dz_required = 0;
    r_plus = NaN;
    r_minus = NaN;
    r_solve = NaN;
    sensitivity = NaN;
    trialCode = uint8(0);

    % rejectCode meaning:
    % 0 = not eligible
    % 1 = already zero residual
    % 2 = invalid finite-difference residuals
    % 3 = sensitivity too small
    % 4 = solved dz invalid or zero
    % 5 = solved trial projection invalid
    % 6 = solved trial did not improve residual
    % 7 = accepted Newton/Jacobi solve
    % 8 = accepted fallback +increment
    % 9 = accepted fallback -increment
    rejectCode = uint8(0);

    if ~pre_eligible || ~isfinite(pre_r0)
        return
    end

    eligible = true;
    p_sel = pre_p_sel;
    opt_uv_idx = pre_opt_uv_idx;
    abs_uv_idx = pre_abs_uv_idx;
    r0 = pre_r0;

    abs0 = abs(r0);
    if ~isfinite(abs0) || abs0 == 0
        trialCode = uint8(1);
        rejectCode = uint8(1);
        return
    end

    % Estimate local residual sensitivity to WSE.
    [r_plus, r_minus] = KLT_test_single_vector_perturbation( ...
        camA, ...
        pre_uvA, ...
        pre_uvB, ...
        dem_current, usePreparedDEMRowTrials, demPrepCurrent, ...
        ix_now, iy_now, increment, ...
        TransX, TransY, phi, shift, flowAngleOffsetDeg);

    % Prefer central difference, fall back to one-sided if needed.
    if isfinite(r_plus) && isfinite(r_minus)
        sensitivity = (r_plus - r_minus) / (2 * increment);
    elseif isfinite(r_plus)
        sensitivity = (r_plus - r0) / increment;
    elseif isfinite(r_minus)
        sensitivity = (r0 - r_minus) / increment;
    else
        sensitivity = NaN;
        rejectCode = uint8(2);
    end

    solvedAccepted = false;

    if isfinite(sensitivity)
        if abs(sensitivity) < jacobiMinSensitivity
            rejectCode = uint8(3);
        else
            dz_try = -r0 / sensitivity;
            dz_try = max(-jacobiMaxStep, min(jacobiMaxStep, dz_try));

            if ~isfinite(dz_try) || dz_try == 0
                rejectCode = uint8(4);
            else
                if jacobiRequireImprovement
                    r_solve = KLT_test_single_signed_perturbation( ...
                        camA, ...
                        pre_uvA, ...
                        pre_uvB, ...
                        dem_current, usePreparedDEMRowTrials, demPrepCurrent, ...
                        ix_now, iy_now, dz_try, ...
                        TransX, TransY, phi, shift, flowAngleOffsetDeg);

                    if ~isfinite(r_solve)
                        rejectCode = uint8(5);
                    elseif abs(r_solve) < abs0
                        dz_required = dz_try;
                        trialCode = uint8(2);
                        rejectCode = uint8(7);
                        solvedAccepted = true;
                    else
                        rejectCode = uint8(6);
                    end
                else
                    dz_required = dz_try;
                    trialCode = uint8(2);
                    rejectCode = uint8(7);
                    solvedAccepted = true;
                end
            end
        end
    end

    if solvedAccepted
        return
    end

    % Robust fallback to the original +/- increment choice.
    if jacobiFallbackToIncrement
        absPlus = Inf;
        absMinus = Inf;

        if isfinite(r_plus)
            absPlus = abs(r_plus);
        end
        if isfinite(r_minus)
            absMinus = abs(r_minus);
        end

        if absPlus < abs0 && absPlus <= absMinus
            dz_required = +increment;
            r_solve = r_plus;
            trialCode = uint8(3);
            rejectCode = uint8(8);
        elseif absMinus < abs0 && absMinus < absPlus
            dz_required = -increment;
            r_solve = r_minus;
            trialCode = uint8(4);
            rejectCode = uint8(9);
        else
            dz_required = 0;
            trialCode = uint8(1);
            if rejectCode == 0
                rejectCode = uint8(6);
            end
        end
    else
        dz_required = 0;
        trialCode = uint8(1);
        if rejectCode == 0
            rejectCode = uint8(6);
        end
    end
end

% =========================================================================
% Convert Jacobi trial code / reject code to log text
% =========================================================================
function txt = KLT_jacobi_trial_code_to_text(code)
    switch uint8(code)
        case 1
            txt = 'none';
        case 2
            txt = 'Jacobi solve';
        case 3
            txt = 'fallback +increment';
        case 4
            txt = 'fallback -increment';
        otherwise
            txt = 'N/A';
    end
end

function txt = KLT_reject_code_to_text(code)
    switch uint8(code)
        case 0
            txt = 'not eligible';
        case 1
            txt = 'already zero';
        case 2
            txt = 'invalid finite difference';
        case 3
            txt = 'low sensitivity';
        case 4
            txt = 'invalid dz';
        case 5
            txt = 'trial projection invalid';
        case 6
            txt = 'trial did not improve';
        case 7
            txt = 'accepted Newton';
        case 8
            txt = 'accepted fallback +';
        case 9
            txt = 'accepted fallback -';
        otherwise
            txt = 'unknown';
    end
end

% =========================================================================
% Test both +-increment perturbations on ONE vector
% =========================================================================
function [r_plus, r_minus] = KLT_test_single_vector_perturbation( ...
    camA, uv_A_single, uv_B_single, ...
    dem_current, usePreparedDEMRowTrials, demPrepCurrent, ...
    ix_now, iy_now, increment, ...
    TransX, TransY, phi, shift, flowAngleOffsetDeg)

    r_plus = NaN;
    r_minus = NaN;

    if usePreparedDEMRowTrials
        prep_plus = KLT_build_single_cell_trial_prep( ...
            camA, demPrepCurrent, dem_current, iy_now, ix_now, +increment);
        xyzA_p = camA.invproject(uv_A_single, prep_plus);
        xyzB_p = camA.invproject(uv_B_single, prep_plus);
    else
        dem_plus = KLT_single_cell_perturbation(dem_current, iy_now, ix_now, +increment);
        xyzA_p = camA.invproject(uv_A_single, TransX, TransY, dem_plus);
        xyzB_p = camA.invproject(uv_B_single, TransX, TransY, dem_plus);
    end

    psi_p = KLT_compute_residual_from_projected( ...
        xyzA_p(:,1:2), xyzB_p(:,1:2), phi, shift, flowAngleOffsetDeg);
    if ~isempty(psi_p)
        r_plus = psi_p(1);
    end

    if usePreparedDEMRowTrials
        prep_minus = KLT_build_single_cell_trial_prep( ...
            camA, demPrepCurrent, dem_current, iy_now, ix_now, -increment);
        xyzA_m = camA.invproject(uv_A_single, prep_minus);
        xyzB_m = camA.invproject(uv_B_single, prep_minus);
    else
        dem_minus = KLT_single_cell_perturbation(dem_current, iy_now, ix_now, -increment);
        xyzA_m = camA.invproject(uv_A_single, TransX, TransY, dem_minus);
        xyzB_m = camA.invproject(uv_B_single, TransX, TransY, dem_minus);
    end

    psi_m = KLT_compute_residual_from_projected( ...
        xyzA_m(:,1:2), xyzB_m(:,1:2), phi, shift, flowAngleOffsetDeg);
    if ~isempty(psi_m)
        r_minus = psi_m(1);
    end
end

% =========================================================================
% Test a SINGLE signed perturbation on ONE vector
% =========================================================================
function r_out = KLT_test_single_signed_perturbation( ...
    camA, uv_A_single, uv_B_single, ...
    dem_current, usePreparedDEMRowTrials, demPrepCurrent, ...
    ix_now, iy_now, dz_test, ...
    TransX, TransY, phi, shift, flowAngleOffsetDeg)

    r_out = NaN;

    if usePreparedDEMRowTrials
        prep_test = KLT_build_single_cell_trial_prep( ...
            camA, demPrepCurrent, dem_current, iy_now, ix_now, dz_test);
        xyzA_t = camA.invproject(uv_A_single, prep_test);
        xyzB_t = camA.invproject(uv_B_single, prep_test);
    else
        dem_test = KLT_single_cell_perturbation(dem_current, iy_now, ix_now, dz_test);
        xyzA_t = camA.invproject(uv_A_single, TransX, TransY, dem_test);
        xyzB_t = camA.invproject(uv_B_single, TransX, TransY, dem_test);
    end

    psi_t = KLT_compute_residual_from_projected( ...
        xyzA_t(:,1:2), xyzB_t(:,1:2), phi, shift, flowAngleOffsetDeg);
    if ~isempty(psi_t)
        r_out = psi_t(1);
    end
end

% =========================================================================
% Advanced convergence diagnostics
% =========================================================================
function metrics = KLT_ensure_advanced_convergence_metrics(metrics, maxIter)

    scalarFields = {
        'endIterAllPathMeanAbsResidualHist'
        'endIterAllPathMedianAbsResidualHist'
        'endIterAllPathRMSResidualHist'
        'endIterAllPathMeanSignedResidualHist'
        'endIterAllPathMedianSignedResidualHist'
        'endIterAllPathFiniteFracHist'

        'endIterRepCellMeanAbsResidualHist'
        'endIterRepCellMedianAbsResidualHist'
        'endIterRepCellRMSResidualHist'
        'endIterRepCellMeanSignedResidualHist'
        'endIterRepCellMedianSignedResidualHist'
        'endIterRepCellCountHist'

        'fixedPathMeanAbsResidualHist'
        'fixedPathMedianAbsResidualHist'
        'fixedPathRMSResidualHist'
        'fixedPathCoverageHist'

        'fixedRepPathMeanAbsResidualHist'
        'fixedRepPathMedianAbsResidualHist'
        'fixedRepPathRMSResidualHist'
        'fixedRepPathCoverageHist'

        'localImprovementFracHist'
        'localWorseningFracHist'
        'localNoChangeFracHist'
        'localMedianImprovementHist'
        'localMeanImprovementHist'
        'localP25ImprovementHist'
        'localP75ImprovementHist'
        'localMedianRelativeImprovementHist'
        'localMeanRelativeImprovementHist'
        'localFracReducedBy10PercentHist'
        'localFracWorseBy10PercentHist'

        'nChangedCellsHist'
        'changedCellFracHist'
        'meanAbsSurfaceChangeHist'
        'medianAbsSurfaceChangeHist'
        'p95AbsSurfaceChangeHist'
        'maxAbsSurfaceChangeHist'
        'netSurfaceChangeHist'
        'positiveUpdateFracHist'
        'negativeUpdateFracHist'

        'signFlipFracHist'
        'sameSignFracHist'

        'activeJaccardHist'
        'activeGainedCellsHist'
        'activeLostCellsHist'

        'fracAcceptedNewtonHist'
        'fracAcceptedFallbackPlusHist'
        'fracAcceptedFallbackMinusHist'
        'fracRejectedLowSensitivityHist'
        'fracRejectedInvalidFDHist'
        'fracRejectedInvalidDzHist'
        'fracRejectedInvalidTrialHist'
        'fracRejectedNoImprovementHist'
        'fracAlreadyZeroHist'

        'convergenceFlagHist'
        'stalledFlagHist'
    };

    for k = 1:numel(scalarFields)
        f = scalarFields{k};

        old = [];
        if isfield(metrics, f) && isnumeric(metrics.(f))
            old = metrics.(f);
        end

        tmp = nan(maxIter, 1);

        if ~isempty(old)
            old = old(:);
            nCopy = min(numel(old), maxIter);
            tmp(1:nCopy) = old(1:nCopy);
        end

        metrics.(f) = tmp;
    end

    matrixFields = {
        'bandMedianAbsResidualHist'
        'bandMeanAbsResidualHist'
        'bandRMSResidualHist'
        'bandUpdatedCellCountHist'
        'bandMedianAbsDzHist'
    };

    for k = 1:numel(matrixFields)
        f = matrixFields{k};
        if ~isfield(metrics, f) || ~isnumeric(metrics.(f))
            metrics.(f) = nan(maxIter, 5);
        else
            old = metrics.(f);
            tmp = nan(maxIter, max(5, size(old,2)));
            nRow = min(size(old,1), maxIter);
            nCol = min(size(old,2), size(tmp,2));
            tmp(1:nRow,1:nCol) = old(1:nRow,1:nCol);
            metrics.(f) = tmp;
        end
    end

    cellFields = {
        'activeMaskHist'
        'rejectCodeMapHist'
        'endIterRepResidualMapHist'
        'rejectCodeHist'
    };

    for k = 1:numel(cellFields)
        f = cellFields{k};

        old = {};
        if isfield(metrics, f) && iscell(metrics.(f))
            old = metrics.(f);
        end

        tmp = cell(maxIter, 1);
        if ~isempty(old)
            old = old(:);
            nCopy = min(numel(old), maxIter);
            tmp(1:nCopy) = old(1:nCopy);
        end
        metrics.(f) = tmp;
    end

    if ~isfield(metrics, 'fixedPathIdx')
        metrics.fixedPathIdx = [];
    end
    if ~isfield(metrics, 'fixedRepPathIdx')
        metrics.fixedRepPathIdx = [];
    end
end

function [repResidualMap, repSelectedLocalIdx, valid_path] = ...
    KLT_build_end_iteration_representative_map( ...
    xyzA_end, xyzB_end, psi_end_all, ...
    xi, yi, dx, dy, origin, e_s, e_n, ...
    Ny, Nx, enforceCurrentDiffCol)

    repResidualMap = nan(Ny, Nx);
    repSelectedLocalIdx = nan(Ny, Nx);
    valid_path = false(size(psi_end_all));

    if isempty(psi_end_all)
        return
    end

    xA = xyzA_end(:,1);
    yA = xyzA_end(:,2);
    xB = xyzB_end(:,1);
    yB = xyzB_end(:,2);

    [cellA_x, cellA_y] = KLT_points_to_cells(xA, yA, xi, yi, dx, dy, origin, e_s, e_n);
    [cellB_x, cellB_y] = KLT_points_to_cells(xB, yB, xi, yi, dx, dy, origin, e_s, e_n);

    valid_bins = isfinite(cellA_x) & isfinite(cellA_y) & ...
        isfinite(cellB_x) & isfinite(cellB_y);

    finite_xy = isfinite(xA) & isfinite(yA) & isfinite(xB) & isfinite(yB);

    valid_geom = isfinite(psi_end_all) & finite_xy & valid_bins;
    valid_diffcol = valid_geom & (cellA_x ~= cellB_x);

    if enforceCurrentDiffCol
        valid_path = valid_diffcol;
    else
        valid_path = valid_geom;
    end

    [~, cell_selected_idx, ~] = KLT_select_one_path_per_cell( ...
        xA, yA, psi_end_all, valid_path, ...
        xi, yi, dx, dy, origin, e_s, e_n);

    repSelectedLocalIdx = cell_selected_idx;
    good = isfinite(cell_selected_idx);

    if any(good(:))
        idx_sel = round(cell_selected_idx(good));
        repResidualMap(good) = psi_end_all(idx_sel);
    end
end

function metrics = KLT_update_advanced_convergence_diagnostics( ...
    metrics, a, maxIter, ...
    residual_used_pre, residual_used_post, ...
    dz_iter, dzMapHist, wse_prev, wse_curr, ...
    nSweepCells, activeMaskThisIter, ...
    y_order, x_order, Ny, Nx, ...
    psi_end_all, endRepResidualMap, endRepSelectedLocalIdx, endValidPathMask, ...
    selectedOptIdxThisIter, reject_code_iter, ...
    nConvergenceBands, convergenceWindow)

    metrics = KLT_ensure_advanced_convergence_metrics(metrics, maxIter);

    % ---------------------------------------------------------------------
    % True end-of-iteration all-path residuals
    % ---------------------------------------------------------------------
    finite_end = isfinite(psi_end_all);
    psi_good = psi_end_all(finite_end);

    if ~isempty(psi_good)
        metrics.endIterAllPathMeanAbsResidualHist(a) = mean(abs(psi_good), 'omitnan');
        metrics.endIterAllPathMedianAbsResidualHist(a) = median(abs(psi_good), 'omitnan');
        metrics.endIterAllPathRMSResidualHist(a) = sqrt(mean(psi_good.^2, 'omitnan'));
        metrics.endIterAllPathMeanSignedResidualHist(a) = mean(psi_good, 'omitnan');
        metrics.endIterAllPathMedianSignedResidualHist(a) = median(psi_good, 'omitnan');
    end
    metrics.endIterAllPathFiniteFracHist(a) = nnz(finite_end) / max(numel(psi_end_all), 1);

    % ---------------------------------------------------------------------
    % True end-of-iteration representative-cell residuals
    % ---------------------------------------------------------------------
    r_rep = endRepResidualMap(isfinite(endRepResidualMap));
    if ~isempty(r_rep)
        metrics.endIterRepCellMeanAbsResidualHist(a) = mean(abs(r_rep), 'omitnan');
        metrics.endIterRepCellMedianAbsResidualHist(a) = median(abs(r_rep), 'omitnan');
        metrics.endIterRepCellRMSResidualHist(a) = sqrt(mean(r_rep.^2, 'omitnan'));
        metrics.endIterRepCellMeanSignedResidualHist(a) = mean(r_rep, 'omitnan');
        metrics.endIterRepCellMedianSignedResidualHist(a) = median(r_rep, 'omitnan');
    end
    metrics.endIterRepCellCountHist(a) = numel(r_rep);
    metrics.endIterRepResidualMapHist{a} = endRepResidualMap;

    % ---------------------------------------------------------------------
    % Fixed-population path residuals
    % ---------------------------------------------------------------------
    if isempty(metrics.fixedPathIdx)
        metrics.fixedPathIdx = find(finite_end);
    end

    if isempty(metrics.fixedRepPathIdx)
        fixedRep = unique(selectedOptIdxThisIter(isfinite(selectedOptIdxThisIter)));
        if isempty(fixedRep)
            fixedRep = unique(endRepSelectedLocalIdx(isfinite(endRepSelectedLocalIdx)));
        end
        metrics.fixedRepPathIdx = fixedRep(:);
    end

    if ~isempty(metrics.fixedPathIdx)
        idxFixed = metrics.fixedPathIdx(:);
        idxFixed = idxFixed(idxFixed >= 1 & idxFixed <= numel(psi_end_all));
        psi_fixed = psi_end_all(idxFixed);
        good_fixed = isfinite(psi_fixed);

        if any(good_fixed)
            metrics.fixedPathMeanAbsResidualHist(a) = mean(abs(psi_fixed(good_fixed)), 'omitnan');
            metrics.fixedPathMedianAbsResidualHist(a) = median(abs(psi_fixed(good_fixed)), 'omitnan');
            metrics.fixedPathRMSResidualHist(a) = sqrt(mean(psi_fixed(good_fixed).^2, 'omitnan'));
        end
        metrics.fixedPathCoverageHist(a) = nnz(good_fixed) / max(numel(idxFixed), 1);
    end

    if ~isempty(metrics.fixedRepPathIdx)
        idxFixedRep = metrics.fixedRepPathIdx(:);
        idxFixedRep = idxFixedRep(idxFixedRep >= 1 & idxFixedRep <= numel(psi_end_all));
        psi_fixed_rep = psi_end_all(idxFixedRep);
        good_fixed_rep = isfinite(psi_fixed_rep);

        if any(good_fixed_rep)
            metrics.fixedRepPathMeanAbsResidualHist(a) = mean(abs(psi_fixed_rep(good_fixed_rep)), 'omitnan');
            metrics.fixedRepPathMedianAbsResidualHist(a) = median(abs(psi_fixed_rep(good_fixed_rep)), 'omitnan');
            metrics.fixedRepPathRMSResidualHist(a) = sqrt(mean(psi_fixed_rep(good_fixed_rep).^2, 'omitnan'));
        end
        metrics.fixedRepPathCoverageHist(a) = nnz(good_fixed_rep) / max(numel(idxFixedRep), 1);
    end

    % ---------------------------------------------------------------------
    % Local improvement distributions
    % ---------------------------------------------------------------------
    if ~isempty(residual_used_pre) && ~isempty(residual_used_post)
        n = min(numel(residual_used_pre), numel(residual_used_post));
        absPre = abs(residual_used_pre(1:n));
        absPost = abs(residual_used_post(1:n));

        localImprovement = absPre - absPost;
        goodImp = isfinite(localImprovement);

        if any(goodImp)
            imp = localImprovement(goodImp);
            metrics.localImprovementFracHist(a) = mean(imp > 0, 'omitnan');
            metrics.localWorseningFracHist(a) = mean(imp < 0, 'omitnan');
            metrics.localNoChangeFracHist(a) = mean(imp == 0, 'omitnan');
            metrics.localMedianImprovementHist(a) = median(imp, 'omitnan');
            metrics.localMeanImprovementHist(a) = mean(imp, 'omitnan');
            metrics.localP25ImprovementHist(a) = prctile(imp, 25);
            metrics.localP75ImprovementHist(a) = prctile(imp, 75);
        end

        goodRel = isfinite(absPre) & isfinite(absPost) & absPre > 0;
        if any(goodRel)
            relativeImprovement = (absPre(goodRel) - absPost(goodRel)) ./ absPre(goodRel);
            metrics.localMedianRelativeImprovementHist(a) = median(relativeImprovement, 'omitnan');
            metrics.localMeanRelativeImprovementHist(a) = mean(relativeImprovement, 'omitnan');
            metrics.localFracReducedBy10PercentHist(a) = mean(relativeImprovement > 0.10, 'omitnan');
            metrics.localFracWorseBy10PercentHist(a) = mean(relativeImprovement < -0.10, 'omitnan');
        end
    end

    % ---------------------------------------------------------------------
    % WSE surface-change diagnostics
    % ---------------------------------------------------------------------
    if ~isempty(wse_prev) && ~isempty(wse_curr) && isequal(size(wse_prev), size(wse_curr))
        dz_surface = wse_curr - wse_prev;
        dz_changed = dz_surface(isfinite(dz_surface) & dz_surface ~= 0);

        metrics.nChangedCellsHist(a) = numel(dz_changed);
        metrics.changedCellFracHist(a) = numel(dz_changed) / max(nSweepCells, 1);

        if ~isempty(dz_changed)
            metrics.meanAbsSurfaceChangeHist(a) = mean(abs(dz_changed), 'omitnan');
            metrics.medianAbsSurfaceChangeHist(a) = median(abs(dz_changed), 'omitnan');
            metrics.p95AbsSurfaceChangeHist(a) = prctile(abs(dz_changed), 95);
            metrics.maxAbsSurfaceChangeHist(a) = max(abs(dz_changed));
            metrics.netSurfaceChangeHist(a) = sum(dz_changed, 'omitnan');
            metrics.positiveUpdateFracHist(a) = mean(dz_changed > 0, 'omitnan');
            metrics.negativeUpdateFracHist(a) = mean(dz_changed < 0, 'omitnan');
        end
    end

    % ---------------------------------------------------------------------
    % Sign oscillation diagnostics
    % ---------------------------------------------------------------------
    if a > 1 && numel(dzMapHist) >= a-1 && ~isempty(dzMapHist{a-1})
        dz_prev = dzMapHist{a-1};
        dz_curr = dz_iter;

        if isequal(size(dz_prev), size(dz_curr))
            common = isfinite(dz_prev) & isfinite(dz_curr) & dz_prev ~= 0 & dz_curr ~= 0;
            if any(common(:))
                metrics.signFlipFracHist(a) = mean(sign(dz_prev(common)) ~= sign(dz_curr(common)), 'omitnan');
                metrics.sameSignFracHist(a) = mean(sign(dz_prev(common)) == sign(dz_curr(common)), 'omitnan');
            end
        end
    end

    % ---------------------------------------------------------------------
    % Active-population churn diagnostics
    % ---------------------------------------------------------------------
    metrics.activeMaskHist{a} = activeMaskThisIter;

    if a > 1 && numel(metrics.activeMaskHist) >= a-1 && ~isempty(metrics.activeMaskHist{a-1})
        activePrev = metrics.activeMaskHist{a-1};
        activeNow = activeMaskThisIter;

        if isequal(size(activePrev), size(activeNow))
            unionMask = activePrev | activeNow;
            interMask = activePrev & activeNow;

            if any(unionMask(:))
                metrics.activeJaccardHist(a) = nnz(interMask) / nnz(unionMask);
            end
            metrics.activeGainedCellsHist(a) = nnz(~activePrev & activeNow);
            metrics.activeLostCellsHist(a) = nnz(activePrev & ~activeNow);
        end
    end

    % ---------------------------------------------------------------------
    % Rejection / acceptance reason diagnostics
    % ---------------------------------------------------------------------
    metrics.rejectCodeMapHist{a} = reject_code_iter;

    rejectCodes = reject_code_iter(isfinite(reject_code_iter));
    metrics.rejectCodeHist{a} = uint8(rejectCodes(:));

    if ~isempty(rejectCodes)
        metrics.fracAcceptedNewtonHist(a) = mean(rejectCodes == 7, 'omitnan');
        metrics.fracAcceptedFallbackPlusHist(a) = mean(rejectCodes == 8, 'omitnan');
        metrics.fracAcceptedFallbackMinusHist(a) = mean(rejectCodes == 9, 'omitnan');
        metrics.fracRejectedLowSensitivityHist(a) = mean(rejectCodes == 3, 'omitnan');
        metrics.fracRejectedInvalidFDHist(a) = mean(rejectCodes == 2, 'omitnan');
        metrics.fracRejectedInvalidDzHist(a) = mean(rejectCodes == 4, 'omitnan');
        metrics.fracRejectedInvalidTrialHist(a) = mean(rejectCodes == 5, 'omitnan');
        metrics.fracRejectedNoImprovementHist(a) = mean(rejectCodes == 6, 'omitnan');
        metrics.fracAlreadyZeroHist(a) = mean(rejectCodes == 1, 'omitnan');
    end

    % ---------------------------------------------------------------------
    % Row-band convergence diagnostics
    % ---------------------------------------------------------------------
    nBands = max(1, nConvergenceBands);
    if size(metrics.bandMedianAbsResidualHist,2) < nBands
        metrics.bandMedianAbsResidualHist(:,end+1:nBands) = NaN;
        metrics.bandMeanAbsResidualHist(:,end+1:nBands) = NaN;
        metrics.bandRMSResidualHist(:,end+1:nBands) = NaN;
        metrics.bandUpdatedCellCountHist(:,end+1:nBands) = NaN;
        metrics.bandMedianAbsDzHist(:,end+1:nBands) = NaN;
    end

    if ~isempty(y_order) && ~isempty(x_order)
        rowMin = min(y_order);
        rowMax = max(y_order);
        edges = round(linspace(rowMin, rowMax + 1, nBands + 1));

        for b = 1:nBands
            rows_b = edges(b):(edges(b+1)-1);
            rows_b = rows_b(rows_b >= 1 & rows_b <= Ny);

            if isempty(rows_b)
                continue
            end

            mask_b = false(Ny, Nx);
            mask_b(rows_b, x_order) = true;

            r_b = endRepResidualMap(mask_b);
            r_b = r_b(isfinite(r_b));

            dz_b = dz_iter(mask_b);
            dz_b = dz_b(isfinite(dz_b) & dz_b ~= 0);

            if ~isempty(r_b)
                metrics.bandMedianAbsResidualHist(a,b) = median(abs(r_b), 'omitnan');
                metrics.bandMeanAbsResidualHist(a,b) = mean(abs(r_b), 'omitnan');
                metrics.bandRMSResidualHist(a,b) = sqrt(mean(r_b.^2, 'omitnan'));
            end

            metrics.bandUpdatedCellCountHist(a,b) = numel(dz_b);

            if ~isempty(dz_b)
                metrics.bandMedianAbsDzHist(a,b) = median(abs(dz_b), 'omitnan');
            end
        end
    end

    % ---------------------------------------------------------------------
    % Diagnostic convergence / stall flags
    % These flags do NOT stop the solver. Thresholds are intentionally exposed.
    % ---------------------------------------------------------------------
    metrics.convergenceThresholds.window = convergenceWindow;
    metrics.convergenceThresholds.residualStableDeg = 0.1;
    metrics.convergenceThresholds.smallMedianDzM = 0.002;
    metrics.convergenceThresholds.fewChangedCellsFrac = 0.01;
    metrics.convergenceThresholds.activeJaccard = 0.95;
    metrics.convergenceThresholds.signFlipFrac = 0.30;
    metrics.convergenceThresholds.lowResidualDeg = 1.0;

    if a >= convergenceWindow
        recent = a-convergenceWindow+1:a;

        residTrend = metrics.fixedRepPathMedianAbsResidualHist(recent);
        if all(~isfinite(residTrend))
            residTrend = metrics.endIterRepCellMedianAbsResidualHist(recent);
        end

        dzTrend = metrics.medianAbsSurfaceChangeHist(recent);
        changeTrend = metrics.changedCellFracHist(recent);

        residStable = KLT_nan_range(residTrend) < metrics.convergenceThresholds.residualStableDeg;
        smallDz = median(dzTrend, 'omitnan') < metrics.convergenceThresholds.smallMedianDzM;
        fewUpdates = median(changeTrend, 'omitnan') < metrics.convergenceThresholds.fewChangedCellsFrac;

        if all(isfinite(metrics.activeJaccardHist(recent)))
            stablePopulation = median(metrics.activeJaccardHist(recent), 'omitnan') > metrics.convergenceThresholds.activeJaccard;
        else
            stablePopulation = true;
        end

        signVals = metrics.signFlipFracHist(recent);
        if any(isfinite(signVals))
            notOscillating = median(signVals, 'omitnan') < metrics.convergenceThresholds.signFlipFrac;
        else
            notOscillating = true;
        end

        lowResidual = median(residTrend, 'omitnan') < metrics.convergenceThresholds.lowResidualDeg;

        metrics.convergenceFlagHist(a) = double(residStable && smallDz && fewUpdates && stablePopulation && notOscillating && lowResidual);
        metrics.stalledFlagHist(a) = double(residStable && smallDz && fewUpdates && ~lowResidual);
    end

    % Useful final snapshots
    metrics.endIterAllPathMedianAbsResidualFinal = KLT_last_finite(metrics.endIterAllPathMedianAbsResidualHist);
    metrics.endIterRepCellMedianAbsResidualFinal = KLT_last_finite(metrics.endIterRepCellMedianAbsResidualHist);
    metrics.fixedRepPathMedianAbsResidualFinal = KLT_last_finite(metrics.fixedRepPathMedianAbsResidualHist);
    metrics.medianAbsSurfaceChangeFinal = KLT_last_finite(metrics.medianAbsSurfaceChangeHist);
    metrics.changedCellFracFinal = KLT_last_finite(metrics.changedCellFracHist);
    metrics.activeJaccardFinal = KLT_last_finite(metrics.activeJaccardHist);
    metrics.signFlipFracFinal = KLT_last_finite(metrics.signFlipFracHist);
end

% =========================================================================
% Standard convergence metric updater
% =========================================================================
function metrics = KLT_update_convergence_metrics( ...
    metrics, ...
    meanAbsResidualHist, medianAbsResidualHist, rmsResidualHist, ...
    meanSignedResidualHist, medianSignedResidualHist, ...
    meanAbsResidualPostHist, medianAbsResidualPostHist, rmsResidualPostHist, ...
    meanSignedResidualPostHist, medianSignedResidualPostHist, ...
    medianAbsDzHist, maxAbsDzHist)

    metrics.meanAbsResidualHist = meanAbsResidualHist;
    metrics.medianAbsResidualHist = medianAbsResidualHist;
    metrics.rmsResidualHist = rmsResidualHist;
    metrics.meanSignedResidualHist = meanSignedResidualHist;
    metrics.medianSignedResidualHist = medianSignedResidualHist;

    metrics.meanAbsResidualPostHist = meanAbsResidualPostHist;
    metrics.medianAbsResidualPostHist = medianAbsResidualPostHist;
    metrics.rmsResidualPostHist = rmsResidualPostHist;
    metrics.meanSignedResidualPostHist = meanSignedResidualPostHist;
    metrics.medianSignedResidualPostHist = medianSignedResidualPostHist;

    metrics.medianAbsDzHist = medianAbsDzHist;
    metrics.maxAbsDzHist = maxAbsDzHist;

    % Rename conceptually in interpretation as "within-sweep" improvement.
    metrics.withinSweepMeanAbsImprovementHist = ...
        meanAbsResidualHist - meanAbsResidualPostHist;

    metrics.withinSweepMedianAbsImprovementHist = ...
        medianAbsResidualHist - medianAbsResidualPostHist;

    metrics.withinSweepRMSImprovementHist = ...
        rmsResidualHist - rmsResidualPostHist;

    metrics.withinSweepMeanAbsResidualRatioHist = ...
        KLT_safe_ratio(meanAbsResidualPostHist, meanAbsResidualHist);

    metrics.withinSweepMedianAbsResidualRatioHist = ...
        KLT_safe_ratio(medianAbsResidualPostHist, medianAbsResidualHist);

    metrics.withinSweepRMSResidualRatioHist = ...
        KLT_safe_ratio(rmsResidualPostHist, rmsResidualHist);

    % Backwards-compatible names
    metrics.meanAbsResidualImprovementHist = metrics.withinSweepMeanAbsImprovementHist;
    metrics.medianAbsResidualImprovementHist = metrics.withinSweepMedianAbsImprovementHist;
    metrics.rmsResidualImprovementHist = metrics.withinSweepRMSImprovementHist;
    metrics.meanAbsResidualRatioHist = metrics.withinSweepMeanAbsResidualRatioHist;
    metrics.medianAbsResidualRatioHist = metrics.withinSweepMedianAbsResidualRatioHist;
    metrics.rmsResidualRatioHist = metrics.withinSweepRMSResidualRatioHist;

    metrics.meanAbsResidualFinal = KLT_last_finite(meanAbsResidualHist);
    metrics.medianAbsResidualFinal = KLT_last_finite(medianAbsResidualHist);
    metrics.rmsResidualFinal = KLT_last_finite(rmsResidualHist);
    metrics.meanSignedResidualFinal = KLT_last_finite(meanSignedResidualHist);
    metrics.medianSignedResidualFinal = KLT_last_finite(medianSignedResidualHist);

    metrics.meanAbsResidualPostFinal = KLT_last_finite(meanAbsResidualPostHist);
    metrics.medianAbsResidualPostFinal = KLT_last_finite(medianAbsResidualPostHist);
    metrics.rmsResidualPostFinal = KLT_last_finite(rmsResidualPostHist);
    metrics.meanSignedResidualPostFinal = KLT_last_finite(meanSignedResidualPostHist);
    metrics.medianSignedResidualPostFinal = KLT_last_finite(medianSignedResidualPostHist);

    metrics.withinSweepMeanAbsImprovementFinal = ...
        KLT_last_finite(metrics.withinSweepMeanAbsImprovementHist);

    metrics.withinSweepMedianAbsImprovementFinal = ...
        KLT_last_finite(metrics.withinSweepMedianAbsImprovementHist);

    metrics.withinSweepRMSImprovementFinal = ...
        KLT_last_finite(metrics.withinSweepRMSImprovementHist);

    metrics.medianAbsDzFinal = KLT_last_finite(medianAbsDzHist);
    metrics.maxAbsDzFinal = KLT_last_finite(maxAbsDzHist);

    metrics.meanAbsResidualAllIters = mean(meanAbsResidualHist, 'omitnan');
    metrics.medianAbsResidualAllIters = mean(medianAbsResidualHist, 'omitnan');
    metrics.rmsResidualAllIters = mean(rmsResidualHist, 'omitnan');
    metrics.meanSignedResidualAllIters = mean(meanSignedResidualHist, 'omitnan');
    metrics.medianSignedResidualAllIters = mean(medianSignedResidualHist, 'omitnan');

    metrics.meanAbsResidualPostAllIters = mean(meanAbsResidualPostHist, 'omitnan');
    metrics.medianAbsResidualPostAllIters = mean(medianAbsResidualPostHist, 'omitnan');
    metrics.rmsResidualPostAllIters = mean(rmsResidualPostHist, 'omitnan');
    metrics.meanSignedResidualPostAllIters = mean(meanSignedResidualPostHist, 'omitnan');
    metrics.medianSignedResidualPostAllIters = mean(medianSignedResidualPostHist, 'omitnan');

    metrics.medianAbsDzAllIters = mean(medianAbsDzHist, 'omitnan');
    metrics.maxAbsDzAllIters = mean(maxAbsDzHist, 'omitnan');
end

function out = KLT_safe_ratio(num, den)
    out = nan(size(num));
    good = isfinite(num) & isfinite(den) & den ~= 0;
    out(good) = num(good) ./ den(good);
end

function r = KLT_nan_range(x)
    x = x(isfinite(x));
    if isempty(x)
        r = Inf;
    else
        r = max(x) - min(x);
    end
end

% =========================================================================
% Infer flow-angle offset from saved/initial residual metrics
% =========================================================================
function [offsetDeg, medianRawDeg] = KLT_infer_flow_angle_offset_from_metrics(metrics, rowRange, colRange)
    offsetDeg = 0;
    medianRawDeg = NaN;

    if nargin < 1 || isempty(metrics) || ~isstruct(metrics)
        return
    end
    if nargin < 2
        rowRange = [];
    end
    if nargin < 3
        colRange = [];
    end

    % Prefer the raw initial residual map when it is available. This is
    % important when resuming from a checkpoint, because initCellMedianVal
    % may already have had the flow-angle offset removed.
    if isfield(metrics, 'initCellMedianValRaw') && ~isempty(metrics.initCellMedianValRaw)
        residualMapRaw = metrics.initCellMedianValRaw;
    elseif isfield(metrics, 'initCellMedianVal') && ~isempty(metrics.initCellMedianVal)
        residualMapRaw = metrics.initCellMedianVal;
    else
        return
    end

    [offsetDeg, medianRawDeg] = ...
        KLT_median_residual_in_window(residualMapRaw, rowRange, colRange);
end

function [offsetDeg, medianRawDeg] = KLT_median_residual_in_window(residualMapRaw, rowRange, colRange)
    offsetDeg = 0;
    medianRawDeg = NaN;

    if isempty(residualMapRaw)
        return
    end
    if nargin < 2
        rowRange = [];
    end
    if nargin < 3
        colRange = [];
    end

    sz = size(residualMapRaw);
    Ny = sz(1);
    if numel(sz) >= 2
        Nx = sz(2);
    else
        Nx = 1;
    end

    if isempty(rowRange)
        rows = 1:Ny;
    else
        rows = round(min(rowRange)):round(max(rowRange));
        rows = rows(rows >= 1 & rows <= Ny);
        if isempty(rows)
            warning('Requested residual zeroing row range does not overlap residual map rows 1:%d.', Ny);
            return
        end
    end

    if isempty(colRange)
        cols = 1:Nx;
    else
        cols = round(min(colRange)):round(max(colRange));
        cols = cols(cols >= 1 & cols <= Nx);
        if isempty(cols)
            warning('Requested residual zeroing column range does not overlap residual map columns 1:%d.', Nx);
            return
        end
    end

    subs = repmat({':'}, 1, ndims(residualMapRaw));
    subs{1} = rows;
    if ndims(residualMapRaw) >= 2
        subs{2} = cols;
    end
    vals = residualMapRaw(subs{:});
    vals = vals(isfinite(vals));

    if isempty(vals)
        if isempty(rowRange) && isempty(colRange)
            warning('No finite residual values found across the initial residual domain.');
        elseif ~isempty(rowRange) && ~isempty(colRange)
            warning('No finite residual values found in requested residual zeroing rows and columns.');
        elseif ~isempty(rowRange)
            warning('No finite residual values found in requested residual zeroing rows.');
        else
            warning('No finite residual values found in requested residual zeroing columns.');
        end
        return
    end

    medianRawDeg = median(vals(:), 'omitnan');
    if isfinite(medianRawDeg)
        offsetDeg = medianRawDeg;
    end
end

% =========================================================================
% Compute angular residual psi from projected path endpoints
% =========================================================================
function psi = KLT_compute_residual_from_projected(xyzA, xyzB, phi, shift, flowAngleOffsetDeg)
    if nargin < 5 || isempty(flowAngleOffsetDeg)
        flowAngleOffsetDeg = 0;
    end

    v2 = xyzB - xyzA;
    obs_dir = rad2deg(atan2(v2(:,2), v2(:,1)));
    obs_dir = KLT_applyAngleShift(obs_dir, shift);

    obs_dir = obs_dir + flowAngleOffsetDeg;
    psi = phi - obs_dir;
end

% =========================================================================
% Angle helpers
% =========================================================================
function [ang_out, shift] = KLT_wrapTo360_centerMedian(ang_in)
    % This simple version preserves the current angle convention.
    % Replace with your original centre-median version if it used different wrapping.
    shift = 0;
    ang_out = mod(ang_in + shift, 360);
end

function ang_out = KLT_applyAngleShift(ang_in, shift)
    if nargin < 2 || isempty(shift)
        shift = 0;
    end
    ang_out = mod(ang_in + shift, 360);
end

% =========================================================================
% Convert invproject output into a single [x y] point
% =========================================================================
function xy = KLT_force_xy_point(p)
    if isempty(p)
        error('invproject returned an empty result when estimating flow direction.');
    end

    sz = size(p);
    if isvector(p)
        p = p(:).';
        if numel(p) < 2
            error('invproject vector output has fewer than 2 elements.');
        end
        xy = p(1:2);
        return
    end

    if ismatrix(p)
        if sz(2) >= 2
            xy = p(1,1:2);
            return
        elseif sz(1) >= 2
            xy = p(1:2,1).';
            return
        end
    end

    if ndims(p) == 3
        if sz(3) >= 2
            xy = [p(1,1,1), p(1,1,2)];
            return
        end
    end

    error('Could not interpret invproject output of size [%s] as an XY point.', num2str(sz));
end

% =========================================================================
% Convert real-world points to rotated-grid cell indices
% =========================================================================
function [ix, iy] = KLT_points_to_cells(x, y, xi, yi, dx, dy, origin, e_s, e_n)
    rel = [x(:) - origin(1), y(:) - origin(2)];
    s = rel * e_s(:);
    n = rel * e_n(:);

    x0 = xi(1) - 0.5*dx;
    y0 = yi(1) - 0.5*dy;

    ix = floor((s - x0) ./ dx) + 1;
    iy = floor((n - y0) ./ dy) + 1;

    in = ix >= 1 & ix <= numel(xi) & iy >= 1 & iy <= numel(yi);
    ix(~in) = NaN;
    iy(~in) = NaN;
end

% =========================================================================
% Select one representative path per cell using A-point assignment
% =========================================================================
function [selected_path_mask, cell_selected_idx, cell_median_val] = ...
    KLT_select_one_path_per_cell( ...
    xA_path, yA_path, path_val, valid_mask, ...
    xi, yi, dx, dy, origin, e_s, e_n)

    nPath = numel(path_val);
    Nx = numel(xi);
    Ny = numel(yi);

    selected_path_mask = false(nPath,1);
    cell_selected_idx = nan(Ny, Nx);
    cell_median_val = nan(Ny, Nx);

    if ~any(valid_mask)
        return
    end

    cell_paths = cell(Ny, Nx);
    cell_vals = cell(Ny, Nx);

    for k = find(valid_mask(:))'
        xa = xA_path(k);
        ya = yA_path(k);
        vp = path_val(k);

        if ~isfinite(xa) || ~isfinite(ya) || ~isfinite(vp)
            continue
        end

        [ix_a, iy_a] = KLT_points_to_cells(xa, ya, xi, yi, dx, dy, origin, e_s, e_n);
        if ~isfinite(ix_a) || ~isfinite(iy_a)
            continue
        end

        lin = sub2ind([Ny, Nx], iy_a, ix_a);
        cell_paths{lin}(end+1) = k; %#ok<AGROW>
        cell_vals{lin}(end+1) = vp; %#ok<AGROW>
    end

    for r = 1:Ny
        for c = 1:Nx
            vals = cell_vals{r,c};
            ids = cell_paths{r,c};
            if isempty(vals)
                continue
            end

            medv = median(vals, 'omitnan');
            if ~isfinite(medv)
                continue
            end

            d = abs(vals - medv);
            d(~isfinite(d)) = inf;
            [bestDist, ii] = min(d); %#ok<ASGLU>
            if isempty(ii) || ~isfinite(bestDist)
                continue
            end

            cell_selected_idx(r,c) = ids(ii);
            cell_median_val(r,c) = medv;
            selected_path_mask(ids(ii)) = true;
        end
    end
end

% =========================================================================
% Select one representative path from a candidate set
% =========================================================================
function [p_sel, medv] = KLT_select_representative_path_for_cell(path_ids, r_path)
    p_sel = NaN;
    medv = NaN;

    if isempty(path_ids)
        return
    end

    vals = r_path(path_ids);
    good = isfinite(vals);
    if ~any(good)
        return
    end

    path_ids = path_ids(good);
    vals = vals(good);
    medv = median(vals, 'omitnan');
    if ~isfinite(medv)
        return
    end

    d = abs(vals - medv);
    d(~isfinite(d)) = inf;
    [bestDist, k] = min(d); %#ok<ASGLU>
    if isempty(k) || ~isfinite(bestDist)
        return
    end

    p_sel = path_ids(k);
end

function out = KLT_last_finite(x)
    idx = find(isfinite(x), 1, 'last');
    if isempty(idx)
        out = NaN;
    else
        out = x(idx);
    end
end

function idx = KLT_last_valid_iter(C)
    idx = find(~cellfun(@isempty, C), 1, 'last');
    if isempty(idx)
        idx = 1;
    end
end

% =========================================================================
% Preselect the representative path for one row on the frozen column state
% =========================================================================
function [eligible, p_sel, opt_uv_idx, abs_uv_idx, r0, uvA, uvB] = ...
    KLT_preselect_row_candidate_scalar( ...
    ids_here, r_path_frozen, abs_uv_for_local, opt_uv_for_local, ...
    xyzA_opt_uv, xyzB_opt_uv)

    eligible = false;
    p_sel = NaN;
    opt_uv_idx = NaN;
    abs_uv_idx = NaN;
    r0 = NaN;
    uvA = nan(1, size(xyzA_opt_uv, 2));
    uvB = nan(1, size(xyzB_opt_uv, 2));

    good_here = isfinite(r_path_frozen(ids_here));
    ids_here = ids_here(good_here);
    if isempty(ids_here)
        return
    end

    [p_sel_tmp, ~] = KLT_select_representative_path_for_cell(ids_here, r_path_frozen);
    if ~isfinite(p_sel_tmp)
        return
    end

    eligible = true;
    p_sel = p_sel_tmp;
    opt_uv_idx = opt_uv_for_local(p_sel_tmp);
    abs_uv_idx = abs_uv_for_local(p_sel_tmp);
    r0 = r_path_frozen(p_sel_tmp);
    uvA = xyzA_opt_uv(opt_uv_idx, :);
    uvB = xyzB_opt_uv(opt_uv_idx, :);
end

% =========================================================================
% DEM perturbation helpers
% =========================================================================
function dem_out = KLT_single_cell_perturbation(dem_current, iy_now, ix_now, dz)
    dem_out = dem_current;
    dem_out(iy_now, ix_now) = dem_out(iy_now, ix_now) + dz;
end

function prep_out = KLT_build_single_cell_trial_prep( ...
    camA, demPrepCurrent, dem_current, iy_now, ix_now, dz)

    dem_trial = KLT_single_cell_perturbation(dem_current, iy_now, ix_now, dz);
    prep_out = camA.updatePreparedDEMInterpOnly(demPrepCurrent, dem_trial);
end

% =========================================================================
% Rebuild A/B cell membership from the current projected cache
% =========================================================================
function [base_idx_local, ...
    abs_uv_for_local, ...
    opt_uv_for_local, ...
    r_path, ...
    cellA_path_ids, ...
    cellB_path_ids, ...
    active_rows_by_col, ...
    nActive] = KLT_rebuild_AB_membership_from_cache( ...
    xyzA_cache, xyzB_cache, psi_cache, opt_idx_full, ...
    xi, yi, dx, dy, origin, e_s, e_n, ...
    Ny, Nx, x_order, y_order, enforceCurrentDiffCol)

    cellA_path_ids = cell(Ny, Nx);
    cellB_path_ids = cell(Ny, Nx);
    active_rows_by_col = cell(1, Nx);
    nActive = 0;

    base_idx_local = find(isfinite(psi_cache));
    abs_uv_for_local = opt_idx_full(base_idx_local);
    opt_uv_for_local = base_idx_local;
    r_path = psi_cache(base_idx_local);

    if isempty(base_idx_local)
        return
    end

    xA_path = xyzA_cache(base_idx_local, 1);
    yA_path = xyzA_cache(base_idx_local, 2);
    xB_path = xyzB_cache(base_idx_local, 1);
    yB_path = xyzB_cache(base_idx_local, 2);

    [cellA_x, cellA_y] = KLT_points_to_cells( ...
        xA_path, yA_path, xi, yi, dx, dy, origin, e_s, e_n);

    [cellB_x, cellB_y] = KLT_points_to_cells( ...
        xB_path, yB_path, xi, yi, dx, dy, origin, e_s, e_n);

    valid_bins = isfinite(cellA_x) & isfinite(cellA_y) & ...
        isfinite(cellB_x) & isfinite(cellB_y);

    valid_geom = isfinite(r_path) & ...
        isfinite(xA_path) & isfinite(yA_path) & ...
        isfinite(xB_path) & isfinite(yB_path) & ...
        valid_bins;

    valid_diffcol = valid_geom & (cellA_x ~= cellB_x);

    if enforceCurrentDiffCol
        valid_path = valid_diffcol;
    else
        valid_path = valid_geom;
    end

    valid_idx = find(valid_path);

    if ~isempty(valid_idx)
        linA_ids = sub2ind([Ny, Nx], cellA_y(valid_idx), cellA_x(valid_idx));
        [linA_sorted, orderA] = sort(linA_ids);
        valid_idx_A_sorted = valid_idx(orderA);
        cutA = [1; find(diff(linA_sorted)) + 1; numel(linA_sorted) + 1];

        for g = 1:numel(cutA)-1
            i1 = cutA(g);
            i2 = cutA(g+1) - 1;
            cellA_path_ids{linA_sorted(i1)} = valid_idx_A_sorted(i1:i2).';
        end

        linB_ids = sub2ind([Ny, Nx], cellB_y(valid_idx), cellB_x(valid_idx));
        [linB_sorted, orderB] = sort(linB_ids);
        valid_idx_B_sorted = valid_idx(orderB);
        cutB = [1; find(diff(linB_sorted)) + 1; numel(linB_sorted) + 1];

        for g = 1:numel(cutB)-1
            i1 = cutB(g);
            i2 = cutB(g+1) - 1;
            cellB_path_ids{linB_sorted(i1)} = valid_idx_B_sorted(i1:i2).';
        end
    end

    for ix_tmp = x_order
        mask_here = ~cellfun('isempty', cellA_path_ids(y_order, ix_tmp));
        rows_here = int32(y_order(mask_here));
        active_rows_by_col{ix_tmp} = rows_here(:);
        nActive = nActive + numel(rows_here);
    end
end

% =========================================================================
% Switchable paired reprojection for A/B UV endpoints
% =========================================================================
function [xyzA, xyzB] = KLT_invproject_uv_pairs_safe_or_prepared( ...
    camA, uvA, uvB, X_rot, Y_rot, dem_current, ...
    usePreparedDEM, demPrepCurrent, xyA0, xyB0)

    nA = size(uvA,1);
    nB = size(uvB,1);

    if nA ~= nB
        error('uvA and uvB must have the same number of rows.');
    end

    if ~usePreparedDEM
        xyzA = camA.invproject(uvA, X_rot, Y_rot, dem_current);
        xyzB = camA.invproject(uvB, X_rot, Y_rot, dem_current);
        return
    end

    uv_pair = [uvA; uvB];
    useWarm = (nargin >= 10) && ~isempty(xyA0) && ~isempty(xyB0);

    if useWarm
        if size(xyA0,1) ~= nA || size(xyB0,1) ~= nB || ...
                size(xyA0,2) ~= 2 || size(xyB0,2) ~= 2
            error('xyA0/xyB0 must be N-by-2 and match uvA/uvB row counts.');
        end

        xy0_pair = [xyA0; xyB0];
        xyz_pair = camA.invproject(uv_pair, demPrepCurrent, xy0_pair);
    else
        xyz_pair = camA.invproject(uv_pair, demPrepCurrent);
    end

    xyzA = xyz_pair(1:nA, :);
    xyzB = xyz_pair(nA+1:end, :);

    bad = any(~isfinite(xyzA(:,1:2)), 2) | any(~isfinite(xyzB(:,1:2)), 2);

    if any(bad)
        xyzA_safe = camA.invproject(uvA(bad, :), X_rot, Y_rot, dem_current);
        xyzB_safe = camA.invproject(uvB(bad, :), X_rot, Y_rot, dem_current);

        xyzA(bad, :) = xyzA_safe;
        xyzB(bad, :) = xyzB_safe;
    end
end

% =========================================================================
% Ensure debug fields exist in metrics
% =========================================================================
function metrics = KLT_ensure_debug_metrics(metrics, maxIter)

    if ~isfield(metrics, 'debug') || ~isstruct(metrics.debug)
        metrics.debug = struct;
    end

    scalarFields = {
        'iterGeomValid'
        'iterDiffColValid'
        'iterLostToSameCol'
        'iterLostToOutOfGrid'
        'iterActiveCells'
        'iterNearActiveCells'
        'iterNearActiveDataCells'
        'iterEligibleCells'
        'iterUpdatedCells'
        'iterTotalOpt'
        'iterFiniteAxy'
        'iterFiniteBxy'
        'iterFiniteBothXY'
        'iterFinitePsi'
        'iterDataRowMin'
        'iterDataRowMax'
        'iterDataSupportedCells'
        'iterLowBandActive'
        'iterHighBandActive'
        'iterLowBandEligible'
        'iterLowBandUpdated'
        'iterHighBandEligible'
        'iterHighBandUpdated'
        'dataRowMinHist'
        'dataRowMaxHist'
        'nDataSupportedCellsHist'
        'activeFracDataSupportedHist'
    };

    for k = 1:numel(scalarFields)
        f = scalarFields{k};

        old = [];
        if isfield(metrics.debug, f) && isnumeric(metrics.debug.(f))
            old = metrics.debug.(f);
        end

        tmp = nan(maxIter, 1);

        if ~isempty(old)
            old = old(:);
            nCopy = min(numel(old), maxIter);
            tmp(1:nCopy) = old(1:nCopy);
        end

        metrics.debug.(f) = tmp;
    end

    cellFields = {
        'candidateCountAStartHist'
        'candidateCountBStartHist'
        'candidateCountAEndHist'
        'candidateCountBEndHist'
        'rowActiveCountStartHist'
        'eligibleNoUpdateHist'
        'nearBandDataHist'
        'dataSupportedMaskHist'
        'activeLowBandHist'
        'activeHighBandHist'
    };

    for k = 1:numel(cellFields)
        f = cellFields{k};

        old = {};
        if isfield(metrics.debug, f) && iscell(metrics.debug.(f))
            old = metrics.debug.(f);
        end

        tmp = cell(maxIter, 1);

        if ~isempty(old)
            old = old(:);
            nCopy = min(numel(old), maxIter);
            tmp(1:nCopy) = old(1:nCopy);
        end

        metrics.debug.(f) = tmp;
    end
end

function [countA, countB] = KLT_count_cell_path_ids(cellA_path_ids, cellB_path_ids)
    [Ny, Nx] = size(cellA_path_ids);
    countA = zeros(Ny, Nx, 'uint16');
    countB = zeros(Ny, Nx, 'uint16');

    for rr = 1:Ny
        for cc = 1:Nx
            countA(rr,cc) = uint16(min(numel(cellA_path_ids{rr,cc}), double(intmax('uint16'))));
            countB(rr,cc) = uint16(min(numel(cellB_path_ids{rr,cc}), double(intmax('uint16'))));
        end
    end
end

function nearActive = KLT_count_near_active_rows(active_rows_by_col, x_order, nearBand)
    nearActive = 0;
    nearBand = double(nearBand(:));

    for ix_tmp = x_order
        rows_here = active_rows_by_col{ix_tmp};
        nearActive = nearActive + nnz(ismember(double(rows_here), nearBand));
    end
end

% =========================================================================
% Load sweep limits from CSV text file using manually supplied video number
% These same effective rows/columns also define the residual-angle calibration window.
% =========================================================================
function [col_start, col_end, row_near, row_far] = ...
    KLT_load_sweep_limits_from_text(configFile, videoNumberForSweep)

    if ~isfile(configFile)
        error('Sweep config file not found: %s', configFile);
    end

    T = readtable(configFile, 'TextType', 'string');
    requiredVars = ["videoNumber", "col_start", "col_end", "row_near", "row_far"];
    missingVars = requiredVars(~ismember(requiredVars, string(T.Properties.VariableNames)));

    if ~isempty(missingVars)
        error('Sweep config file is missing required columns: %s', ...
            strjoin(cellstr(missingVars), ', '));
    end

    videoCol = T.videoNumber;
    if iscell(videoCol)
        videoCol = string(videoCol);
    end

    if isstring(videoCol) || ischar(videoCol)
        match = str2double(string(videoCol)) == double(videoNumberForSweep);
    else
        match = double(videoCol) == double(videoNumberForSweep);
    end

    if ~any(match)
        error('No sweep limits found in %s for videoNumber = %d', ...
            configFile, videoNumberForSweep);
    end

    if nnz(match) > 1
        error('Multiple sweep limit rows found in %s for videoNumber = %d', ...
            configFile, videoNumberForSweep);
    end

    row = T(match, :);
    col_start = double(row.col_start);
    col_end = double(row.col_end);
    row_near = double(row.row_near);
    row_far = double(row.row_far);
end
