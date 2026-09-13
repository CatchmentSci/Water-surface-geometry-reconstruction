function results = KLT_synthetic_truth_test_v3(inputMatFile, testConfig)
% =========================================================================
% Synthetic-truth test harness for KLT_wse_solver_paths_Apoint_block_switchtest_v2
%
% Constructs a known truth WSE field (streamwise sine wave), forward-projects
% synthetic straight streamwise trajectories through it onto the camera,
% then runs the solver from a flat initial DEM and compares the recovered
% surface against truth.
%
% Inputs:
%   inputMatFile - path to the same MAT file that drives the real solver.
%                  Must contain: app_in, camA_fullmodel, camA_first_fullmodel,
%                  xyzA_wse, xyzB_wse, aa, wse_map.
%
%   testConfig   - optional struct with fields:
%                  .caseCsvFile   (char) optional CSV path. When supplied,
%                                 each CSV row is one test case and must contain
%                                 amplitude, wavelength, and either wseRowRange
%                                 or wseRowStart/wseRowEnd columns.
%                  .caseCsvRow    (scalar) optional 1-based CSV data-row index.
%                                 Set this to execute only one CSV row. An empty
%                                 value executes every CSV row (default).
%                  .amplitudes    (vector, m) legacy non-CSV mode only
%                  .wavelengths   (vector, m) legacy non-CSV mode only
%                  .nPathsTarget  (scalar)     default: 100000
%                  .ab_step       (scalar, m)  default: 0.5
%                  .saveDir       (char)       default: pwd
%                  .verbose       (logical)    default: true
%                  .resultsTag    (char)       default: 'syn'
%
% Output:
%   results - struct array, one element per executed CSV row or legacy case.
% =========================================================================

% -------------------------------------------------------------------------
% Defaults
% -------------------------------------------------------------------------
if nargin < 2 || isempty(testConfig)
    testConfig = struct;
elseif ischar(testConfig) || (isstring(testConfig) && isscalar(testConfig))
    % Convenience form: KLT_synthetic_truth_test_v3(inputMatFile, 'cases.csv')
    testConfig = struct('caseCsvFile', char(testConfig));
elseif ~isstruct(testConfig)
    error('testConfig must be a struct or a CSV filename.');
end

% CSV case source and row selection.
% Prefer values supplied by the Slurm wrapper. This makes the CSV row
% consistent with KLT_S_ID, while still allowing explicit testConfig fields
% to override the environment for local/manual runs.
defaultCaseCsvFile = '/mnt/nfs/home/nmp65/Downloads/v1_current/Syn/synthetic_cases.csv';

if ~isfield(testConfig, 'caseCsvFile') || isempty(testConfig.caseCsvFile)
    envCaseCsvFile = strtrim(getenv('KLT_CASE_CSV_FILE'));
    if ~isempty(envCaseCsvFile)
        testConfig.caseCsvFile = envCaseCsvFile;
    else
        testConfig.caseCsvFile = defaultCaseCsvFile;
    end
end

if ~isfield(testConfig, 'caseCsvRow') || isempty(testConfig.caseCsvRow)
    envCaseCsvRow = str2double(getenv('KLT_CASE_CSV_ROW'));
    if isfinite(envCaseCsvRow) && envCaseCsvRow > 0
        testConfig.caseCsvRow = round(envCaseCsvRow);
    else
        envSId = str2double(getenv('KLT_S_ID'));
        if isfinite(envSId) && envSId > 0
            testConfig.caseCsvRow = round(envSId);
        else
            testConfig.caseCsvRow = [];
        end
    end
end

if ~isfield(testConfig, 'caseCsvFile'),  testConfig.caseCsvFile  = '';                 end
if ~isfield(testConfig, 'caseCsvRow'),   testConfig.caseCsvRow   = [];                 end
if ~isfield(testConfig, 'amplitudes'),    testConfig.amplitudes    = [];               end %[0.02, 0.06, 0.10]
if ~isfield(testConfig, 'wavelengths'),   testConfig.wavelengths   = [];               end % replaced from residual-angle 2D FFT by default
if ~isfield(testConfig, 'useResidualFFTForWavelengths'), testConfig.useResidualFFTForWavelengths = false; end
if ~isfield(testConfig, 'nPathsTarget'),  testConfig.nPathsTarget  = 2000000;          end
if ~isfield(testConfig, 'ab_step'),       testConfig.ab_step       = 0.05;             end % insensitive
if ~isfield(testConfig, 'saveDir'),       testConfig.saveDir       = pwd;              end
if ~isfield(testConfig, 'verbose'),       testConfig.verbose       = true;             end
if ~isfield(testConfig, 'resultsTag'),    testConfig.resultsTag    = 'syn';            end
if ~isfield(testConfig, 'wseRowRange'),   testConfig.wseRowRange   = [];               end
if ~isfield(testConfig, 'autoFlowAngleOffset'), testConfig.autoFlowAngleOffset = true; end
if ~isfield(testConfig, 'flowAngleOffsetDeg'), testConfig.flowAngleOffsetDeg = [];     end

% Resume / restart controls.
% resumeFromCheckpoint=true means an existing per-case checkpoint is reused
% instead of being deleted. The solver itself restores x_order/y_order from
% the checkpoint, so the resumed run uses the same rows/columns as the
% original checkpoint.
if ~isfield(testConfig, 'resumeFromCheckpoint'), testConfig.resumeFromCheckpoint = true; end % or flase if no resume
if ~isfield(testConfig, 'resumeCaseIdx'),        testConfig.resumeCaseIdx = [];           end % empty means all matching cases
if ~isfield(testConfig, 'resumeCheckpointFile'), testConfig.resumeCheckpointFile = '';              end
if ~isfield(testConfig, 'resumeCaseStateFile'),  testConfig.resumeCaseStateFile = '';              end
if ~isfield(testConfig, 'forceFresh'),           testConfig.forceFresh = false;           end
if ~isfield(testConfig, 'saveCaseState'),        testConfig.saveCaseState = true;         end

amplitudes   = testConfig.amplitudes(:);
wavelengths  = testConfig.wavelengths(:);
caseCsvFile  = char(testConfig.caseCsvFile);
caseCsvRow   = testConfig.caseCsvRow;
useCsvCases  = ~isempty(strtrim(caseCsvFile));
useResidualFFTForWavelengths = testConfig.useResidualFFTForWavelengths;
nPathsTarget = testConfig.nPathsTarget;
ab_step      = testConfig.ab_step;
saveDir      = testConfig.saveDir;
verbose      = testConfig.verbose;
resultsTag   = testConfig.resultsTag;
wseRowRange  = testConfig.wseRowRange;
resumeFromCheckpoint = testConfig.resumeFromCheckpoint;
resumeCaseIdx        = testConfig.resumeCaseIdx;
manualResumeCheckpointFile = char(testConfig.resumeCheckpointFile);
resumeCaseStateFile  = char(testConfig.resumeCaseStateFile);
forceFresh           = testConfig.forceFresh;
saveCaseState        = testConfig.saveCaseState;

if useCsvCases
    allCsvCases = local_read_case_csv(caseCsvFile);
    nCsvRowsAvailable = numel(allCsvCases);

    if isempty(caseCsvRow)
        testCases = allCsvCases;
    else
        if ~isscalar(caseCsvRow) || ~isnumeric(caseCsvRow) || ...
                ~isfinite(caseCsvRow) || caseCsvRow ~= round(caseCsvRow)
            error('testConfig.caseCsvRow must be an empty value or one finite integer scalar.');
        end
        caseCsvRow = round(caseCsvRow);
        if caseCsvRow < 1 || caseCsvRow > nCsvRowsAvailable
            error('Requested testConfig.caseCsvRow=%d, but CSV contains %d data rows.', ...
                caseCsvRow, nCsvRowsAvailable);
        end
        testCases = allCsvCases(caseCsvRow);
    end

    % A CSV-provided wavelength is explicit and must not be replaced by the
    % residual-angle FFT estimate.
    useResidualFFTForWavelengths = false;
    testConfig.useResidualFFTForWavelengths = false;
    testConfig.caseCsvRow = caseCsvRow;
    testConfig.csvRowsAvailable = nCsvRowsAvailable;
else
    testCases = struct([]);
    nCsvRowsAvailable = 0;
end

if ~isfolder(saveDir)
    mkdir(saveDir);
end

% -------------------------------------------------------------------------
% Load the real input MAT
% -------------------------------------------------------------------------
if ~isfile(inputMatFile)
    error('Input MAT file not found: %s', inputMatFile);
end

% Use the input MAT basename in generated output filenames so that cases
% run from different hydraulic-stage inputs cannot overwrite one another.
[~, inputMatBaseName, ~] = fileparts(inputMatFile);
inputMatToken = local_filename_text(inputMatBaseName);
testConfig.inputMatFile = inputMatFile;
testConfig.inputMatToken = inputMatToken;

S = load(inputMatFile);

requiredVars = {'app_in', 'camA_fullmodel', 'camA_first_fullmodel', ...
    'xyzA_wse', 'xyzB_wse', 'aa', 'wse_map'};

for k = 1:numel(requiredVars)
    if ~isfield(S, requiredVars{k})
        error('Input MAT file is missing required variable: %s', requiredVars{k});
    end
end

app = S.app_in;
app.camA       = camera(S.camA_fullmodel);
app.camA_first = camera(S.camA_first_fullmodel);
aa             = S.aa;
wse_map_real   = S.wse_map;

if isfield(S, 'globalPolarity')
    globalPolarity = S.globalPolarity;
else
    globalPolarity = [];
end

% -------------------------------------------------------------------------
% Recover the rotated-grid geometry the solver will build internally.
% This must match exactly what the solver does at lines ~297–340 of
% KLT_wse_solver_paths_Apoint_block_switchtest_v2.
% -------------------------------------------------------------------------
[X_rot, Y_rot, xi, yi, e_s, e_n, origin, hgt_use, dx, dy, phi, shift] = ...
    local_build_rotated_grid(app);

Nx = numel(xi);
Ny = numel(yi);

% Blank or omitted CSV row ranges mean: use the full available WSE grid.
if useCsvCases
    for k = 1:numel(testCases)
        if isempty(testCases(k).wseRowRange) || all(isnan(testCases(k).wseRowRange))
            testCases(k).wseRowRange = [1 Ny];
        end
    end
end

% -------------------------------------------------------------------------
% Estimate initial real-input diagnostics. In legacy mode this may replace
% the wavelength vector. In CSV mode each row has its own WSE row range, so
% diagnostics are recalculated inside the case loop and CSV wavelengths are
% never overwritten.
% -------------------------------------------------------------------------
flatDEM_for_residual = zeros(size(X_rot)) + hgt_use;

if ~useCsvCases
    inputResidualDiagnostics = local_estimate_initial_residual_diagnostics( ...
        app, S.xyzA_wse{aa}, S.xyzB_wse{aa}, ...
        X_rot, Y_rot, flatDEM_for_residual, ...
        xi, yi, dx, dy, origin, e_s, e_n, phi, shift, wseRowRange);

    if useResidualFFTForWavelengths
        if isfinite(inputResidualDiagnostics.fft.peakWavelength2D) && ...
                inputResidualDiagnostics.fft.peakWavelength2D > 0
            wavelengths = inputResidualDiagnostics.fft.peakWavelength2D;
            testConfig.wavelengths = wavelengths;
        elseif isempty(wavelengths)
            warning(['2D FFT residual wavelength estimate failed; ', ...
                     'falling back to 2.5 m.']);
            wavelengths = 2.5;
            testConfig.wavelengths = wavelengths;
        end
    elseif isempty(wavelengths)
        wavelengths = 2.5;
        testConfig.wavelengths = wavelengths;
    end

    testConfig.inputResidualMedianOffsetDeg = inputResidualDiagnostics.medianResidualDeg;
    testConfig.inputResidualFFTpeakWavelength2D = inputResidualDiagnostics.fft.peakWavelength2D;
    testConfig.inputResidualFFTpeakKx = inputResidualDiagnostics.fft.peakKx;
    testConfig.inputResidualFFTpeakKy = inputResidualDiagnostics.fft.peakKy;

    testCases = local_make_legacy_cases(amplitudes, wavelengths, wseRowRange);
end

nCases = numel(testCases);
testConfig.caseDefinitions = testCases;
testConfig.resumeCheckpointFiles = cell(nCases, 1);
testConfig.resumeCaseStateFiles = cell(nCases, 1);

% A selected single CSV row always has a compact value-based filename. If
% all CSV rows are executed together and a pair is repeated, append the CSV
% row only for those repeated pairs so that checkpoint files remain unique.
caseNeedsCsvRowSuffix = false(nCases, 1);
if useCsvCases && nCases > 1
    pairKeys = strings(nCases, 1);
    for k = 1:nCases
        pairKeys(k) = "amp" + string(local_filename_number(testCases(k).amplitude)) + ...
            "_w" + string(local_filename_number(testCases(k).wavelength));
    end
    [~, ~, pairGroup] = unique(pairKeys);
    pairCounts = accumarray(pairGroup, 1);
    caseNeedsCsvRowSuffix = pairCounts(pairGroup) > 1;
end
testConfig.caseNeedsCsvRowSuffix = caseNeedsCsvRowSuffix;

if verbose
    fprintf('\n========================================\n');
    fprintf('Synthetic-truth test harness\n');
    fprintf('========================================\n');
    fprintf('Grid: Ny=%d, Nx=%d, dx=%.3f m, dy=%.3f m\n', Ny, Nx, dx, dy);
    fprintf('Base WSE level: %.3f m\n', hgt_use);
    fprintf('AB streamwise step: %.3f m\n', ab_step);
    fprintf('Target synthetic vectors: %d\n', nPathsTarget);
    if useCsvCases
        fprintf('CSV case file: %s\n', caseCsvFile);
        fprintf('CSV data rows available: %d\n', nCsvRowsAvailable);
        if isempty(caseCsvRow)
            fprintf('CSV rows selected: all (%d cases)\n', nCases);
        else
            fprintf('CSV row selected: %d (executing one case)\n', caseCsvRow);
        end
    else
        fprintf('Legacy test cases: %d amplitudes x wavelengths = %d\n', ...
            numel(amplitudes), numel(wavelengths), nCases);
        fprintf('Input residual median over rows %d:%d = %+0.6f deg\n', ...
            min(wseRowRange), max(wseRowRange), inputResidualDiagnostics.medianResidualDeg);
        fprintf('Input residual 2D FFT peak wavelength = %.4f m', ...
            inputResidualDiagnostics.fft.peakWavelength2D);
        fprintf(' (kx=%+.5f cyc/m, ky=%+.5f cyc/m)\n', ...
            inputResidualDiagnostics.fft.peakKx, inputResidualDiagnostics.fft.peakKy);
    end
end

% -------------------------------------------------------------------------
% Loop over test cases. CSV mode executes one case per row; legacy mode uses
% the original amplitude/wavelength Cartesian product.
% -------------------------------------------------------------------------
results = struct([]);

for caseIdx = 1:nCases
    amp = testCases(caseIdx).amplitude;
    wln = testCases(caseIdx).wavelength;
    wseRowRange = testCases(caseIdx).wseRowRange;
    sourceCsvRow = testCases(caseIdx).sourceCsvRow;

    % Keep the active CSV row visible in testConfig while this case runs.
    testConfig.amplitude = amp;
    testConfig.wavelength = wln;
    testConfig.wavelegth = wln; % accepted typo alias for existing callers/CSV headers
    testConfig.wseRowRange = wseRowRange;
    testConfig.activeCaseCsvRow = sourceCsvRow;

    if useCsvCases
        inputResidualDiagnostics = local_estimate_initial_residual_diagnostics( ...
            app, S.xyzA_wse{aa}, S.xyzB_wse{aa}, ...
            X_rot, Y_rot, flatDEM_for_residual, ...
            xi, yi, dx, dy, origin, e_s, e_n, phi, shift, wseRowRange);
    end

        if verbose
            fprintf('\n----------------------------------------\n');
            if useCsvCases
                fprintf('Case %d (CSV row %d): amp=%.15g m, wavelength=%.15g m, WSE rows=%d:%d\n', ...
                    caseIdx, sourceCsvRow, amp, wln, min(wseRowRange), max(wseRowRange));
            else
                fprintf('Case %d: amp=%.15g m, wavelength=%.15g m, WSE rows=%d:%d\n', ...
                    caseIdx, amp, wln, min(wseRowRange), max(wseRowRange));
            end
            fprintf('----------------------------------------\n');
        end

        % -----------------------------------------------------------------
        % Per-case checkpoint and case-state file names.
        %
        % The checkpoint is the solver restart file.
        % The case-state file stores the synthetic UV pairs, truth DEM, and
        % offset used to create this case, so a resumed run cannot accidentally
        % regenerate a different synthetic problem.
        % -----------------------------------------------------------------
        if useCsvCases && caseNeedsCsvRowSuffix(caseIdx)
            caseOutputStem = sprintf('%s_%s_case_amp%s_w%s_csvrow%d', ...
                resultsTag, inputMatToken, local_filename_number(amp), ...
                local_filename_number(wln), sourceCsvRow);
        else
            caseOutputStem = sprintf('%s_%s_case_amp%s_w%s', ...
                resultsTag, inputMatToken, local_filename_number(amp), ...
                local_filename_number(wln));
        end

        defaultCkpt = fullfile(saveDir, [caseOutputStem '_ckpt.mat']);

        ckpt = defaultCkpt;

        % CSV mode always uses the deterministic name above. In legacy mode,
        % retain support for a manually supplied resume checkpoint file.
        if ~useCsvCases && resumeFromCheckpoint && ~isempty(manualResumeCheckpointFile)
            if isempty(resumeCaseIdx) || ismember(caseIdx, resumeCaseIdx)
                ckpt = manualResumeCheckpointFile;
            end
        end

        resumeIdxSelected = isempty(resumeCaseIdx) || ismember(caseIdx, resumeCaseIdx) || ...
            (useCsvCases && ismember(sourceCsvRow, resumeCaseIdx));

        if ~isempty(resumeCaseStateFile) && resumeIdxSelected
            caseStateFile = resumeCaseStateFile;
        else
            caseStateFile = regexprep(ckpt, '_ckpt\.mat$', '_case_state.mat');
        end

        % Store the generated filename on testConfig as requested. For a CSV
        % with multiple rows, resumeCheckpointFiles preserves every filename;
        % resumeCheckpointFile is the currently executing case.
        testConfig.resumeCheckpointFile = ckpt;
        testConfig.resumeCheckpointFiles{caseIdx} = ckpt;
        testConfig.resumeCaseStateFile = caseStateFile;
        testConfig.resumeCaseStateFiles{caseIdx} = caseStateFile;

        resumeSelected = resumeIdxSelected;
        if ~useCsvCases && ~isempty(manualResumeCheckpointFile) && ...
                strcmp(ckpt, manualResumeCheckpointFile)
            resumeSelected = true;
        end
        resumeThisCase = resumeFromCheckpoint && isfile(ckpt) && resumeSelected;

        if forceFresh && isfile(ckpt) && ~resumeThisCase
            delete(ckpt);
        elseif forceFresh && isfile(ckpt) && resumeThisCase
            warning('forceFresh=true was ignored for resumed case %d: %s', caseIdx, ckpt);
        end

        % -----------------------------------------------------------------
        % Build or reload the synthetic case state.
        % -----------------------------------------------------------------
        if resumeThisCase && isfile(caseStateFile)
            Ccase = load(caseStateFile, 'caseState');
            caseState = Ccase.caseState;

            truthDEM = caseState.truthDEM;
            xyzA_wse_syn = caseState.xyzA_wse_syn;
            xyzB_wse_syn = caseState.xyzB_wse_syn;
            wse_map_syn = caseState.wse_map_syn;
            flowAngleOffsetDeg = caseState.flowAngleOffsetDeg;
            caseResidualDiagnostics = caseState.caseResidualDiagnostics;
            n_generated = caseState.n_generated;

            if isfield(caseState, 'inputResidualDiagnostics')
                inputResidualDiagnostics = caseState.inputResidualDiagnostics;
            end

            if isfield(caseState, 'wseRowRange')
                wseRowRange = caseState.wseRowRange;
            end

            if verbose
                fprintf('Resuming case %d from checkpoint:\n  %s\n', caseIdx, ckpt);
                fprintf('Loaded saved synthetic case state:\n  %s\n', caseStateFile);
                fprintf('Saved synthetic row range: %d:%d\n', min(wseRowRange), max(wseRowRange));
            end

        else
            if resumeThisCase && ~isfile(caseStateFile)
                warning(['Resuming from checkpoint but no saved case-state file was found:\n  %s\n', ...
                         'The synthetic case will be regenerated from the current testConfig. ', ...
                         'This is safe only if amplitude, wavelength, nPathsTarget, ab_step, ', ...
                         'and wseRowRange are unchanged.'], caseStateFile);
            end

            % -------------------------------------------------------------
            % Build the truth surface: streamwise sine, flat in cross-stream.
            %
            % z(s,n) = hgt_use + amp * sin(2*pi*s/wavelength)
            % -------------------------------------------------------------
            truthDEM = local_build_streamwise_sine(xi, yi, hgt_use, amp, wln);

            % -------------------------------------------------------------
            % Generate synthetic UV pairs.
            % -------------------------------------------------------------
            [uvA_syn, uvB_syn, n_generated] = local_generate_synthetic_UV_pairs( ...
                app.camA, xi, yi, e_s, e_n, origin, truthDEM, ab_step, nPathsTarget, wseRowRange);

            %comparisons_horizontal_transects.m

            if verbose
                fprintf('Generated %d synthetic UV pairs\n', n_generated);
            end

            % -------------------------------------------------------------
            % Pack into the xyzA_wse / xyzB_wse cell array structure expected
            % by the solver. Replace cell aa with the synthetic pairs.
            % -------------------------------------------------------------
            xyzA_wse_syn     = S.xyzA_wse;
            xyzB_wse_syn     = S.xyzB_wse;
            xyzA_wse_syn{aa} = uvA_syn;
            xyzB_wse_syn{aa} = uvB_syn;

            % Initial DEM is flat at hgt_use, like the real solver does.
            wse_map_syn       = wse_map_real;
            wse_map_syn{aa,1} = zeros(size(X_rot)) + hgt_use;

            % -------------------------------------------------------------
            % Estimate and pass the residual-angle offset for this synthetic case.
            % -------------------------------------------------------------
            caseResidualDiagnostics = local_estimate_initial_residual_diagnostics( ...
                app, uvA_syn, uvB_syn, ...
                X_rot, Y_rot, flatDEM_for_residual, ...
                xi, yi, dx, dy, origin, e_s, e_n, phi, shift, wseRowRange);

            if testConfig.autoFlowAngleOffset
                flowAngleOffsetDeg = caseResidualDiagnostics.medianResidualDeg;
            else
                flowAngleOffsetDeg = testConfig.flowAngleOffsetDeg;
            end

            if isempty(flowAngleOffsetDeg) || ~isfinite(flowAngleOffsetDeg)
                flowAngleOffsetDeg = [];
            end

            if saveCaseState
                caseState = struct;
                caseState.caseIdx = caseIdx;
                caseState.sourceCsvRow = sourceCsvRow;
                caseState.inputMatFile = inputMatFile;
                caseState.inputMatToken = inputMatToken;
                caseState.outputStem = caseOutputStem;
                caseState.amplitude = amp;
                caseState.wavelength = wln;
                caseState.truthDEM = truthDEM;
                caseState.xyzA_wse_syn = xyzA_wse_syn;
                caseState.xyzB_wse_syn = xyzB_wse_syn;
                caseState.wse_map_syn = wse_map_syn;
                caseState.globalPolarity = globalPolarity;
                caseState.flowAngleOffsetDeg = flowAngleOffsetDeg;
                caseState.caseResidualDiagnostics = caseResidualDiagnostics;
                caseState.inputResidualDiagnostics = inputResidualDiagnostics;
                caseState.n_generated = n_generated;
                caseState.wseRowRange = wseRowRange;
                caseState.nPathsTarget = nPathsTarget;
                caseState.ab_step = ab_step;
                caseState.resultsTag = resultsTag;
                save(caseStateFile, 'caseState', '-v7.3');
            end
        end

        if verbose
            fprintf('Case raw residual median over rows %d:%d = %+0.6f deg; applying offset %+0.6f deg\n', ...
                min(wseRowRange), max(wseRowRange), ...
                caseResidualDiagnostics.medianResidualDeg, ...
                local_printable_scalar(flowAngleOffsetDeg));

            if resumeThisCase
                Sck = load(ckpt, 'checkpoint');
                fprintf('Checkpoint says: lastCompletedIter=%d, nextIter=%d\n', ...
                    Sck.checkpoint.lastCompletedIter, Sck.checkpoint.nextIter);
                fprintf('Checkpoint sweep rows: %d:%d; sweep cols: %d:%d\n', ...
                    min(Sck.checkpoint.y_order), max(Sck.checkpoint.y_order), ...
                    min(Sck.checkpoint.x_order), max(Sck.checkpoint.x_order));
            end
        end

        % -----------------------------------------------------------------
        % Run or resume the solver on the synthetic data.
        %
        % On resume, the solver restores x_order/y_order/active_window_mask
        % from the checkpoint, so the same rows/columns are processed as in
        % the checkpoint, regardless of current testConfig changes.
        % -----------------------------------------------------------------

        [wse_map_out, metrics_out] = KLT_wse_solver_paths_Apoint_block_jacobi_v5( ...
            app, xyzA_wse_syn, xyzB_wse_syn, aa, wse_map_syn, globalPolarity, ...
            ckpt, resumeThisCase, [], ...
            [], [], flowAngleOffsetDeg, ...
            wseRowRange);

        % -----------------------------------------------------------------
        % Compare recovered DEM to truth.
        % -----------------------------------------------------------------
        recoveredDEM = wse_map_out{aa, end};

        if isfield(metrics_out, 'sweepMask')
            evalMask = isfinite(metrics_out.sweepMask);
        else
            evalMask = true(size(truthDEM));
        end

        % Restrict evaluation to the same WSE-map rows where synthetic vectors exist.
        wseRowMask = false(size(truthDEM));

        rowMinEval = min(wseRowRange);
        rowMaxEval = max(wseRowRange);

        rowMinEval = max(1, rowMinEval);
        rowMaxEval = min(size(truthDEM, 1), rowMaxEval);

        wseRowMask(rowMinEval:rowMaxEval, :) = true;

        evalMask = evalMask & wseRowMask;

        errorMap = recoveredDEM - truthDEM;
        errorMap(~evalMask) = NaN;

        rmse = sqrt(mean(errorMap(evalMask).^2, 'omitnan'));
        bias = mean(errorMap(evalMask), 'omitnan');

        % -----------------------------------------------------------------
        % Near-camera vs far-camera bands.
        % -----------------------------------------------------------------
        [rowMin, rowMax] = local_active_row_range(evalMask);

        if isfinite(rowMin) && isfinite(rowMax) && rowMax > rowMin
            mid = round((rowMin + rowMax) / 2);

            nearMask = false(size(evalMask));
            farMask  = false(size(evalMask));

            nearMask(rowMin:mid, :)     = evalMask(rowMin:mid, :);
            farMask(mid+1:rowMax, :)    = evalMask(mid+1:rowMax, :);

            rmse_near = sqrt(mean(errorMap(nearMask).^2, 'omitnan'));
            rmse_far  = sqrt(mean(errorMap(farMask).^2,  'omitnan'));
        else
            rmse_near = NaN;
            rmse_far  = NaN;
        end

        % -----------------------------------------------------------------
        % Streamwise profiles.
        % -----------------------------------------------------------------
        col_mean_err   = nan(1, Nx);
        col_mean_truth = nan(1, Nx);
        col_mean_recov = nan(1, Nx);

        for c = 1:Nx
            col = errorMap(:, c);
            col_mean_err(c) = mean(col, 'omitnan');

            col_t = truthDEM(:, c);
            col_t(~evalMask(:, c)) = NaN;
            col_mean_truth(c) = mean(col_t, 'omitnan');

            col_r = recoveredDEM(:, c);
            col_r(~evalMask(:, c)) = NaN;
            col_mean_recov(c) = mean(col_r, 'omitnan');
        end

        if verbose
            fprintf('RMSE      = %.4f m\n', rmse);
            fprintf('Bias      = %+.4f m\n', bias);
            fprintf('RMSE near = %.4f m\n', rmse_near);
            fprintf('RMSE far  = %.4f m\n', rmse_far);
        end

        % -----------------------------------------------------------------
        % Store results.
        % -----------------------------------------------------------------
        results(caseIdx).sourceCsvRow             = sourceCsvRow; %#ok<AGROW>
        results(caseIdx).inputMatFile             = inputMatFile;
        results(caseIdx).inputMatToken            = inputMatToken;
        results(caseIdx).outputStem               = caseOutputStem;
        results(caseIdx).amplitude               = amp;
        results(caseIdx).wavelength              = wln;
        results(caseIdx).wseRowRange             = wseRowRange;
        results(caseIdx).truthDEM                = truthDEM;
        results(caseIdx).recoveredDEM            = recoveredDEM;
        results(caseIdx).errorMap                = errorMap;
        results(caseIdx).rmse                    = rmse;
        results(caseIdx).bias                    = bias;
        results(caseIdx).rmse_near               = rmse_near;
        results(caseIdx).rmse_far                = rmse_far;
        results(caseIdx).col_mean_err            = col_mean_err;
        results(caseIdx).col_mean_truth          = col_mean_truth;
        results(caseIdx).col_mean_recov          = col_mean_recov;
        results(caseIdx).xi                      = xi;
        results(caseIdx).yi                      = yi;
        results(caseIdx).evalMask                = evalMask;
        results(caseIdx).nSyntheticPaths         = n_generated;
        results(caseIdx).inputResidualDiagnostics = inputResidualDiagnostics;
        results(caseIdx).caseResidualDiagnostics = caseResidualDiagnostics;
        results(caseIdx).flowAngleOffsetDeg      = flowAngleOffsetDeg;
        results(caseIdx).checkpointFile          = ckpt;
        results(caseIdx).caseStateFile           = caseStateFile;
        results(caseIdx).resumedFromCheckpoint   = resumeThisCase;

        if isfield(metrics_out, 'meanAbsResidualHist')
            results(caseIdx).meanAbsResidualHist = metrics_out.meanAbsResidualHist;
        else
            results(caseIdx).meanAbsResidualHist = [];
        end

        if isfield(metrics_out, 'meanAbsResidualPostHist')
            results(caseIdx).meanAbsResidualPostHist = metrics_out.meanAbsResidualPostHist;
        else
            results(caseIdx).meanAbsResidualPostHist = [];
        end

        if isfield(metrics_out, 'activeCellFracHist')
            results(caseIdx).activeCellFracHist = metrics_out.activeCellFracHist;
        else
            results(caseIdx).activeCellFracHist = [];
        end

        if isfield(metrics_out, 'updatedCellFracHist')
            results(caseIdx).updatedCellFracHist = metrics_out.updatedCellFracHist;
        else
            results(caseIdx).updatedCellFracHist = [];
        end

        results(caseIdx).metrics = metrics_out;
end

% -------------------------------------------------------------------------
% Save the full results struct
% -------------------------------------------------------------------------
if useCsvCases && ~isempty(caseCsvRow)
    resultsFile = fullfile(saveDir, sprintf('%s_%s_case_amp%s_w%s_results.mat', ...
        resultsTag, inputMatToken, ...
        local_filename_number(results(1).amplitude), ...
        local_filename_number(results(1).wavelength)));
else
    resultsFile = fullfile(saveDir, sprintf('%s_%s_results.mat', ...
        resultsTag, inputMatToken));
end
save(resultsFile, 'results', 'testConfig', '-v7.3');

if verbose
    fprintf('\n========================================\n');
    fprintf('Test complete. Results saved to:\n  %s\n', resultsFile);
    fprintf('========================================\n');
    fprintf('\nSummary table:\n');
    if useCsvCases
        fprintf('%-4s %-6s %-9s %-10s %-9s %-9s %-9s %-9s %-10s\n', ...
            'case', 'csvRow', 'amp(m)', 'wlen(m)', 'rmse', 'bias', 'rmseNear', 'rmseFar', 'nPaths');
        for k = 1:numel(results)
            fprintf('%-4d %-6d %-9.3f %-10.3f %-9.4f %+9.4f %-9.4f %-9.4f %-10d\n', ...
                k, results(k).sourceCsvRow, results(k).amplitude, results(k).wavelength, ...
                results(k).rmse, results(k).bias, ...
                results(k).rmse_near, results(k).rmse_far, ...
                results(k).nSyntheticPaths);
        end
    else
        fprintf('%-4s %-9s %-10s %-9s %-9s %-9s %-9s %-10s\n', ...
            'case', 'amp(m)', 'wlen(m)', 'rmse', 'bias', 'rmseNear', 'rmseFar', 'nPaths');
        for k = 1:numel(results)
            fprintf('%-4d %-9.3f %-10.3f %-9.4f %+9.4f %-9.4f %-9.4f %-10d\n', ...
                k, results(k).amplitude, results(k).wavelength, ...
                results(k).rmse, results(k).bias, ...
                results(k).rmse_near, results(k).rmse_far, ...
                results(k).nSyntheticPaths);
        end
    end
end

end


% =========================================================================
% Read CSV-defined synthetic cases.
%
% Accepted headers (case-insensitive; punctuation is ignored):
%   amplitude,wavelength,wseRowStart,wseRowEnd
% or
%   amplitude,wavelength,wseRowRange
%
% Example wseRowRange cell values: "210 270", "210:270", or "[210,270]".
% Each CSV row is one case. The main function can select a single CSV row
% with testConfig.caseCsvRow. Duplicate amplitude/wavelength pairs are rejected
% because they would intentionally map to the same requested checkpoint name.
% =========================================================================
function testCases = local_read_case_csv(csvFile)

    if ~isfile(csvFile)
        error('CSV case file not found: %s', csvFile);
    end

    try
        T = readtable(csvFile, 'TextType', 'string', 'VariableNamingRule', 'preserve');
    catch
        % Compatibility fallback for MATLAB releases without VariableNamingRule.
        T = readtable(csvFile);
    end

    if height(T) == 0
        error('CSV case file contains no data rows: %s', csvFile);
    end

    names = T.Properties.VariableNames;
    normalisedNames = cellfun(@local_normalise_csv_header, names, 'UniformOutput', false);

    ampIdx = local_find_csv_column(normalisedNames, ...
        {'amplitude', 'amp', 'testconfigamplitude', 'testconfigamplitudes'}, true, 'amplitude');
    wlnIdx = local_find_csv_column(normalisedNames, ...
        {'wavelength', 'wavelegth', 'wln', 'testconfigwavelength', 'testconfigwavelegth'}, true, 'wavelength');

    rowStartIdx = local_find_csv_column(normalisedNames, ...
        {'wserowstart', 'rowstart', 'wserangestart'}, false, 'wseRowStart');
    rowEndIdx = local_find_csv_column(normalisedNames, ...
        {'wserowend', 'rowend', 'wserangeend'}, false, 'wseRowEnd');
    rowRangeIdx = local_find_csv_column(normalisedNames, ...
        {'wserowrange', 'testconfigwserowrange'}, false, 'wseRowRange');

    amplitudes = local_numeric_table_column(T, ampIdx, 'amplitude');
    wavelengths = local_numeric_table_column(T, wlnIdx, 'wavelength');

    if ~isempty(rowStartIdx) && ~isempty(rowEndIdx)
        rowStart = local_numeric_table_column(T, rowStartIdx, 'wseRowStart');
        rowEnd = local_numeric_table_column(T, rowEndIdx, 'wseRowEnd');
        ranges = [rowStart, rowEnd];
    elseif ~isempty(rowRangeIdx)
        ranges = local_parse_wse_row_ranges(T{:, rowRangeIdx});
    else
        ranges = nan(height(T), 2);
    end

    if any(~isfinite(amplitudes))
        error('CSV amplitude values must be finite.');
    end
    if any(~isfinite(wavelengths) | wavelengths <= 0)
        error('CSV wavelength values must be positive and finite.');
    end

    missingRange = all(isnan(ranges), 2);
    partialMissingRange = any(isnan(ranges), 2) & ~missingRange;
    if any(partialMissingRange)
        error('Each CSV WSE row range must contain either two row indices or be left blank.');
    end

    explicitRanges = ranges(~missingRange, :);
    if any(~isfinite(explicitRanges(:)))
        error('CSV WSE row ranges must be finite when supplied.');
    end
    if any(abs(explicitRanges(:) - round(explicitRanges(:))) > 1e-9)
        error('CSV WSE row ranges must use integer row indices.');
    end

    explicitRanges = round(explicitRanges);
    if any(explicitRanges(:, 1) < 1) || any(explicitRanges(:, 2) <= explicitRanges(:, 1))
        error('Each supplied CSV WSE row range must satisfy 1 <= start < end.');
    end
    ranges(~missingRange, :) = explicitRanges;

    nCases = height(T);
    testCases = repmat(struct('sourceCsvRow', [], 'amplitude', [], 'wavelength', [], 'wseRowRange', []), nCases, 1);

    for k = 1:nCases
        testCases(k).sourceCsvRow = k;
        testCases(k).amplitude = amplitudes(k);
        testCases(k).wavelength = wavelengths(k);
        testCases(k).wseRowRange = ranges(k, :);
    end
end

% =========================================================================
% Construct the original Cartesian-product cases when no CSV is supplied.
% =========================================================================
function testCases = local_make_legacy_cases(amplitudes, wavelengths, wseRowRange)

    nCases = numel(amplitudes) * numel(wavelengths);
    testCases = repmat(struct('sourceCsvRow', [], 'amplitude', [], 'wavelength', [], 'wseRowRange', []), nCases, 1);

    caseIdx = 0;
    for ia = 1:numel(amplitudes)
        for iw = 1:numel(wavelengths)
            caseIdx = caseIdx + 1;
            testCases(caseIdx).sourceCsvRow = NaN;
            testCases(caseIdx).amplitude = amplitudes(ia);
            testCases(caseIdx).wavelength = wavelengths(iw);
            testCases(caseIdx).wseRowRange = wseRowRange;
        end
    end
end

% =========================================================================
% Helpers for robust CSV parsing and deterministic checkpoint filenames.
% =========================================================================
function normalised = local_normalise_csv_header(name)
    normalised = lower(regexprep(char(name), '[^a-zA-Z0-9]', ''));
end

function idx = local_find_csv_column(normalisedNames, candidates, required, displayName)
    idx = [];
    for k = 1:numel(candidates)
        hit = find(strcmp(normalisedNames, candidates{k}), 1, 'first');
        if ~isempty(hit)
            idx = hit;
            return
        end
    end

    if required
        error('CSV is missing required column: %s', displayName);
    end
end

function values = local_numeric_table_column(T, idx, displayName)
    raw = T{:, idx};
    if isnumeric(raw) || islogical(raw)
        values = double(raw);
    else
        values = str2double(string(raw));
    end

    values = values(:);
    if numel(values) ~= height(T) || any(~isfinite(values))
        error('CSV column %s must contain one finite numeric value per row.', displayName);
    end
end

function ranges = local_parse_wse_row_ranges(raw)
    n = numel(raw);
    ranges = nan(n, 2);

    for k = 1:n
        if isnumeric(raw)
            if ~isfinite(raw(k))
                continue
            end
            txt = num2str(raw(k));
        else
            value = string(raw(k));
            if ismissing(value) || strlength(strtrim(value)) == 0
                continue
            end
            txt = char(value);
        end

        % Treat a hyphen between positive digits as a separator, while still
        % allowing signed values to be detected and rejected later.
        txt = regexprep(txt, '(?<=[0-9])\s*-\s*(?=[0-9])', ' ');
        tokens = regexp(txt, '[-+]?(?:\d*\.?\d+)(?:[eE][-+]?\d+)?', 'match');

        if numel(tokens) ~= 2
            error(['Each wseRowRange CSV entry must contain exactly two row ', ...
                   'indices, for example "210 270" or "[210,270]".']);
        end

        ranges(k, :) = [str2double(tokens{1}), str2double(tokens{2})];
    end
end

function token = local_filename_number(x)
    if ~isscalar(x) || ~isfinite(x)
        error('Checkpoint filename values must be finite scalars.');
    end
    token = sprintf('%.15g', x);
end

function token = local_filename_text(txt)
    token = regexprep(char(string(txt)), '[^a-zA-Z0-9_-]+', '_');
    token = regexprep(token, '_+', '_');
    token = regexprep(token, '^_+|_+$', '');
    if isempty(token)
        token = 'input';
    end
end

% =========================================================================
% Build the rotated-grid geometry that the solver builds internally.
% This mirrors lines 297-340 of KLT_wse_solver_paths_Apoint_block_switchtest_v2.
% =========================================================================
function [X_rot, Y_rot, xi, yi, e_s, e_n, origin, hgt_use, dx, dy, phi, shift] = ...
    local_build_rotated_grid(app)

    xIn = app.pts(1,:)';
    yIn = app.pts(2,:)';

    app.start1 = [xIn(1), yIn(1)];
    app.end1   = [xIn(2), yIn(2)];

    frameSize = size(app.firstFrame);
    t1 = frameSize(1);
    t2 = frameSize(2);

    app.start1(1) = max(0, min(t2, app.start1(1)));
    app.start1(2) = max(0, min(t1, app.start1(2)));
    app.end1(1)   = max(0, min(t2, app.end1(1)));
    app.end1(2)   = max(0, min(t1, app.end1(2)));

    hgt_use = app.riverLevelAnalysis(app.videoNumber);

    if ~strcmp(app.OrientationDropDown.Value, 'Dynamic: Stabilisation') && ...
            ~strcmp(app.OrientationDropDown.Value, 'Planet [beta]')

        tempDEM = zeros(size(app.TransX)) + hgt_use;

        start1_rw_raw = app.camA_first.invproject( ...
            app.start1, app.TransX, app.TransY, tempDEM);

        end1_rw_raw = app.camA_first.invproject( ...
            app.end1, app.TransX, app.TransY, tempDEM);

        start1_rw = local_force_xy_point(start1_rw_raw);
        end1_rw   = local_force_xy_point(end1_rw_raw);

    else
        start1_rw = app.start1 .* app.imageResolution;
        end1_rw   = app.end1   .* app.imageResolution;

        start1_rw = start1_rw(1:2);
        end1_rw   = end1_rw(1:2);
    end

    ideal = end1_rw - start1_rw;
    ideal = ideal(:).';

    L_ideal = hypot(ideal(1), ideal(2));

    if L_ideal == 0 || ~isfinite(L_ideal)
        error('Ideal downstream direction is invalid.');
    end

    e_s = ideal / L_ideal;
    e_n = [-e_s(2), e_s(1)];

    try
        [phi, shift] = KLT_wrapTo360_centerMedian( ...
            rad2deg(atan2(ideal(2), ideal(1))));
    catch
        [phi, shift] = local_wrapTo360_centerMedian( ...
            rad2deg(atan2(ideal(2), ideal(1))));
    end

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

    dx = target_dx;
    dy = target_dy;
end


% =========================================================================
% Build a streamwise sine-wave truth surface.
%
% z(s,n) = hgt_use + amp * sin(2*pi*s/wavelength)
% =========================================================================
function truthDEM = local_build_streamwise_sine(xi, yi, hgt_use, amp, wln)

    if wln <= 0 || ~isfinite(wln)
        error('Wavelength must be positive and finite.');
    end

    [Sg, ~] = meshgrid(xi, yi);

    truthDEM = hgt_use + amp * sin(2*pi*Sg/wln);
end


% =========================================================================
% Generate synthetic UV pairs.
%
% Creates approximately nPathsTarget streamwise A->B vectors, but only within
% selected rows of the WSE/checkpoint grid.
%
% wseRowRange = [200 420] means:
%
%   use only yi(200:420)
%
% because rows of checkpoint.wse_map correspond to rows of truthDEM, X_rot,
% Y_rot, and entries in yi.
%
% A and B lie on the truth WSE surface:
%
%   A = (s, n, z_truth(s,n))
%   B = (s + ab_step, n, z_truth(s + ab_step,n))
% =========================================================================
function [uvA, uvB, nGenerated] = local_generate_synthetic_UV_pairs( ...
    camA, xi, yi, e_s, e_n, origin, truthDEM, ab_step, nPathsTarget, wseRowRange)

    if nargin < 9 || isempty(nPathsTarget)
        nPathsTarget = 100000;
    end

    if nargin < 10 || isempty(wseRowRange)
        wseRowRange = [200 420];
    end

    if nPathsTarget <= 0 || ~isfinite(nPathsTarget)
        error('nPathsTarget must be positive and finite.');
    end

    if ab_step <= 0 || ~isfinite(ab_step)
        error('ab_step must be positive and finite.');
    end

    nPathsTarget = round(nPathsTarget);

    % ---------------------------------------------------------------------
    % Convert requested WSE-map row range into yi coordinate limits.
    % ---------------------------------------------------------------------
    Ny = numel(yi);

    rowMinKeep = min(wseRowRange);
    rowMaxKeep = max(wseRowRange);

    rowMinKeep = round(rowMinKeep);
    rowMaxKeep = round(rowMaxKeep);

    if rowMinKeep < 1 || rowMaxKeep > Ny
        error('Requested wseRowRange [%d %d] is outside available WSE rows 1:%d.', ...
            rowMinKeep, rowMaxKeep, Ny);
    end

    if rowMaxKeep <= rowMinKeep
        error('wseRowRange must span at least two rows.');
    end

    % Domain limits in rotated-grid coordinates.
    sMin = min(xi);
    sMax = max(xi) - ab_step;

    % This is the key change:
    % restrict cross-stream coordinate n to the requested WSE-map rows.
    nMin = yi(rowMinKeep);
    nMax = yi(rowMaxKeep);

    if sMax <= sMin
        error('ab_step is too large for the xi domain.');
    end

    sRange = sMax - sMin;
    nRange = nMax - nMin;

    if nRange <= 0
        error('Invalid selected WSE row range: yi(rowMax) must exceed yi(rowMin).');
    end

    % ---------------------------------------------------------------------
    % Build a dense, evenly distributed candidate set only inside rows
    % rowMinKeep:rowMaxKeep of the WSE map.
    % ---------------------------------------------------------------------
    oversampleFactor = 3;
    nCandidateTarget = max(nPathsTarget * oversampleFactor, nPathsTarget);

    % Choose approximately square spacing in physical coordinates.
    nS = ceil(sqrt(nCandidateTarget * sRange / nRange));
    nN = ceil(nCandidateTarget / nS);

    nS = max(nS, 2);
    nN = max(nN, 2);

    sVals = linspace(sMin, sMax, nS);
    nVals = linspace(nMin, nMax, nN);

    [Sg_A, Ng_A] = meshgrid(sVals, nVals);

    % ---------------------------------------------------------------------
    % Deterministic sub-grid staggering.
    % Keeps coverage uniform while avoiding exact alignment with grid nodes.
    % ---------------------------------------------------------------------
    ds = sVals(2) - sVals(1);
    dn = nVals(2) - nVals(1);

    Sg_A = Sg_A + 0.25 * ds * sin(2*pi*(Ng_A - nMin) / max(nRange, eps));
    Ng_A = Ng_A + 0.25 * dn * sin(2*pi*(Sg_A - sMin) / max(sRange, eps));

    % Clip after staggering so points remain inside rows 200:420.
    Sg_A = min(max(Sg_A, sMin), sMax);
    Ng_A = min(max(Ng_A, nMin), nMax);

    % B points are one streamwise step downstream.
    Sg_B = Sg_A + ab_step;
    Ng_B = Ng_A;

    % ---------------------------------------------------------------------
    % Sample truth surface at A and B.
    % ---------------------------------------------------------------------
    Z_interp = griddedInterpolant({yi, xi}, truthDEM, 'linear', 'none');

    Zg_A = Z_interp(Ng_A, Sg_A);
    Zg_B = Z_interp(Ng_B, Sg_B);

    valid = isfinite(Sg_A) & isfinite(Ng_A) & isfinite(Zg_A) & ...
            isfinite(Sg_B) & isfinite(Ng_B) & isfinite(Zg_B);

    Sg_A = Sg_A(valid);
    Ng_A = Ng_A(valid);
    Zg_A = Zg_A(valid);

    Sg_B = Sg_B(valid);
    Ng_B = Ng_B(valid);
    Zg_B = Zg_B(valid);

    % ---------------------------------------------------------------------
    % Convert rotated coordinates back to world XY.
    % ---------------------------------------------------------------------
    Xg_A = origin(1) + Sg_A*e_s(1) + Ng_A*e_n(1);
    Yg_A = origin(2) + Sg_A*e_s(2) + Ng_A*e_n(2);

    Xg_B = origin(1) + Sg_B*e_s(1) + Ng_B*e_n(1);
    Yg_B = origin(2) + Sg_B*e_s(2) + Ng_B*e_n(2);

    xyzA_world = [Xg_A(:), Yg_A(:), Zg_A(:)];
    xyzB_world = [Xg_B(:), Yg_B(:), Zg_B(:)];

    % ---------------------------------------------------------------------
    % Forward-project both endpoints to image coordinates.
    % ---------------------------------------------------------------------
    [uvA_full, ~, inframe_A] = camA.project(xyzA_world);
    [uvB_full, ~, inframe_B] = camA.project(xyzB_world);

    inboth = inframe_A & inframe_B & ...
             all(isfinite(uvA_full), 2) & all(isfinite(uvB_full), 2);

    uvA_all = uvA_full(inboth, :);
    uvB_all = uvB_full(inboth, :);

    nVisible = size(uvA_all, 1);

    if nVisible == 0
        warning(['No synthetic UV pairs survived camera projection. ', ...
            'Camera may not see selected WSE rows %d:%d.'], ...
            rowMinKeep, rowMaxKeep);

        uvA = uvA_all;
        uvB = uvB_all;
        nGenerated = 0;
        return
    end

    % ---------------------------------------------------------------------
    % Keep exactly nPathsTarget if possible.
    % ---------------------------------------------------------------------
    if nVisible >= nPathsTarget
        keepIdx = round(linspace(1, nVisible, nPathsTarget));
        keepIdx = unique(keepIdx, 'stable');

        if numel(keepIdx) < nPathsTarget
            missing = nPathsTarget - numel(keepIdx);
            extraPool = setdiff(1:nVisible, keepIdx, 'stable');
            keepIdx = [keepIdx(:); extraPool(1:missing).'];
        end

        uvA = uvA_all(keepIdx, :);
        uvB = uvB_all(keepIdx, :);
    else
        warning(['Only %d synthetic vectors survived camera projection; ', ...
                 'requested %d. Using all visible vectors.'], ...
                 nVisible, nPathsTarget);

        uvA = uvA_all;
        uvB = uvB_all;
    end

    nGenerated = size(uvA, 1);
end


% =========================================================================
% Mirror of solver's KLT_force_xy_point so the harness is self-contained.
% =========================================================================
function xy = local_force_xy_point(p)

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

    if ndims(p) == 3 && sz(3) >= 2
        xy = [p(1,1,1), p(1,1,2)];
        return
    end

    error('Could not interpret invproject output of size [%s] as an XY point.', ...
        num2str(sz));
end



% =========================================================================
% Estimate initial residual-angle diagnostics from a set of UV path pairs.
% =========================================================================
function diagnostics = local_estimate_initial_residual_diagnostics( ...
    app, uvA, uvB, X_rot, Y_rot, dem_current, ...
    xi, yi, dx, dy, origin, e_s, e_n, phi, shift, wseRowRange)

    diagnostics = struct;
    diagnostics.medianResidualDeg = NaN;
    diagnostics.residualMapRaw = nan(numel(yi), numel(xi));
    diagnostics.residualMapCorrected = nan(numel(yi), numel(xi));
    diagnostics.fft = local_empty_fft_info();

    if isempty(uvA) || isempty(uvB)
        return
    end

    xyzA0 = app.camA.invproject(uvA, X_rot, Y_rot, dem_current);
    xyzB0 = app.camA.invproject(uvB, X_rot, Y_rot, dem_current);

    vel0 = xyzB0 - xyzA0;
    dt = 1;
    try
        if isfinite(app.iter) && isfinite(app.videoFrameRate) && app.videoFrameRate ~= 0
            dt = app.iter * 1 / app.videoFrameRate;
        end
    catch
        dt = 1;
    end

    dvel0 = vel0 ./ dt;
    vmag0 = hypot(dvel0(:,1), dvel0(:,2));
    keep_speed = vmag0 > 0.1;
    idx_speed = find(keep_speed);

    if isempty(idx_speed)
        return
    end

    v20 = xyzB0(idx_speed,1:2) - xyzA0(idx_speed,1:2);
    raw_obs_dir0 = rad2deg(atan2(v20(:,2), v20(:,1)));
    try
        obs_dir0 = KLT_applyAngleShift(raw_obs_dir0, shift);
    catch
        obs_dir0 = local_applyAngleShift(raw_obs_dir0, shift);
    end
    psi0 = phi - obs_dir0;

    filterAngle0 = 60;
    keep_angle = psi0 >= -filterAngle0 & psi0 <= filterAngle0;
    idx_sa = idx_speed(keep_angle);
    psi_sa = psi0(keep_angle);

    if isempty(idx_sa)
        return
    end

    [cellA_x0, cellA_y0] = local_points_to_cells( ...
        xyzA0(idx_sa,1), xyzA0(idx_sa,2), ...
        xi, yi, dx, dy, origin, e_s, e_n);
    [cellB_x0, cellB_y0] = local_points_to_cells( ...
        xyzB0(idx_sa,1), xyzB0(idx_sa,2), ...
        xi, yi, dx, dy, origin, e_s, e_n);

    valid_bins0 = isfinite(cellA_x0) & isfinite(cellA_y0) & ...
        isfinite(cellB_x0) & isfinite(cellB_y0);
    keep_diffcol0 = valid_bins0 & (cellA_x0 ~= cellB_x0);

    idx_cand = idx_sa(keep_diffcol0);
    r_cand = psi_sa(keep_diffcol0);
    valid_mask = true(numel(r_cand), 1);

    if isempty(idx_cand)
        return
    end

    [~, ~, residualMapRaw] = local_select_one_path_per_cell( ...
        xyzA0(idx_cand,1), xyzA0(idx_cand,2), ...
        r_cand, valid_mask, xi, yi, dx, dy, origin, e_s, e_n);

    diagnostics.residualMapRaw = residualMapRaw;
    diagnostics.medianResidualDeg = local_median_residual_in_row_range(residualMapRaw, wseRowRange);

    if isfinite(diagnostics.medianResidualDeg)
        diagnostics.residualMapCorrected = residualMapRaw - diagnostics.medianResidualDeg;
    else
        diagnostics.residualMapCorrected = residualMapRaw;
    end

    diagnostics.fft = local_estimate_2d_peak_wavelength( ...
        diagnostics.residualMapCorrected, xi, yi, wseRowRange);
end

% =========================================================================
% Strongest 2D FFT bin wavelength, not a radial/azimuthal spectral peak.
% =========================================================================
function fftInfo = local_estimate_2d_peak_wavelength(residualMap, xi, yi, rowRange)

    fftInfo = local_empty_fft_info();

    if isempty(residualMap) || isempty(xi) || isempty(yi)
        return
    end

    [Ny, Nx] = size(residualMap);
    if nargin < 4 || isempty(rowRange)
        rows = 1:Ny;
    else
        rows = round(min(rowRange)):round(max(rowRange));
        rows = rows(rows >= 1 & rows <= Ny);
    end

    if numel(rows) < 4 || Nx < 4
        return
    end

    Z = residualMap(rows, :);
    finiteMask = isfinite(Z);

    if nnz(finiteMask) < 4
        return
    end

    globalMedian = median(Z(finiteMask), 'omitnan');
    if ~isfinite(globalMedian)
        globalMedian = 0;
    end

    Zfill = Z;
    for r = 1:size(Zfill,1)
        rowVals = Zfill(r, isfinite(Zfill(r,:)));
        if isempty(rowVals)
            rowFill = globalMedian;
        else
            rowFill = median(rowVals, 'omitnan');
        end
        bad = ~isfinite(Zfill(r,:));
        Zfill(r,bad) = rowFill;
    end

    Zfill = Zfill - mean(Zfill(:), 'omitnan');

    nr = size(Zfill, 1);
    nc = size(Zfill, 2);
    wy = local_cosine_window(nr);
    wx = local_cosine_window(nc).';
    Zwin = Zfill .* (wy * wx);

    P = abs(fftshift(fft2(Zwin))).^2;

    dx_use = median(abs(diff(xi)), 'omitnan');
    dy_use = median(abs(diff(yi(rows))), 'omitnan');

    if ~isfinite(dx_use) || dx_use <= 0 || ~isfinite(dy_use) || dy_use <= 0
        return
    end

    fx = (-floor(nc/2):ceil(nc/2)-1) ./ (nc * dx_use);
    fy = (-floor(nr/2):ceil(nr/2)-1) ./ (nr * dy_use);
    [FX, FY] = meshgrid(fx, fy);
    K = hypot(FX, FY);

    validPeak = isfinite(P) & isfinite(K) & K > 0;
    validIdx = find(validPeak);
    if isempty(validIdx)
        return
    end

    [peakPower, ii] = max(P(validIdx));
    peakLin = validIdx(ii);
    [peakRow, peakCol] = ind2sub(size(P), peakLin);

    peakKx = FX(peakLin);
    peakKy = FY(peakLin);
    peakK = K(peakLin);

    if ~isfinite(peakK) || peakK <= 0
        return
    end

    fftInfo.peakWavelength2D = 1 ./ peakK;
    fftInfo.peakKx = peakKx;
    fftInfo.peakKy = peakKy;
    fftInfo.peakPower = peakPower;
    fftInfo.peakRow = peakRow;
    fftInfo.peakCol = peakCol;
    fftInfo.dx = dx_use;
    fftInfo.dy = dy_use;
    fftInfo.nRows = nr;
    fftInfo.nCols = nc;
    fftInfo.method = 'single strongest 2D FFT bin after DC removal; no radial averaging';

    if peakKx ~= 0
        fftInfo.peakWavelengthX = 1 ./ abs(peakKx);
    else
        fftInfo.peakWavelengthX = Inf;
    end

    if peakKy ~= 0
        fftInfo.peakWavelengthY = 1 ./ abs(peakKy);
    else
        fftInfo.peakWavelengthY = Inf;
    end
end

% =========================================================================
% Empty FFT info struct.
% =========================================================================
function fftInfo = local_empty_fft_info()
    fftInfo = struct( ...
        'peakWavelength2D', NaN, ...
        'peakWavelengthX', NaN, ...
        'peakWavelengthY', NaN, ...
        'peakKx', NaN, ...
        'peakKy', NaN, ...
        'peakPower', NaN, ...
        'peakRow', NaN, ...
        'peakCol', NaN, ...
        'dx', NaN, ...
        'dy', NaN, ...
        'nRows', NaN, ...
        'nCols', NaN, ...
        'method', '');
end

% =========================================================================
% Hann-style cosine taper without toolbox dependencies.
% =========================================================================
function w = local_cosine_window(n)
    if n <= 1
        w = 1;
    else
        w = 0.5 - 0.5*cos(2*pi*(0:n-1)'/(n-1));
    end
end

% =========================================================================
% Median residual in requested WSE rows.
% =========================================================================
function medv = local_median_residual_in_row_range(residualMap, rowRange)
    medv = NaN;
    if isempty(residualMap)
        return
    end

    [Ny, ~] = size(residualMap);
    if nargin < 2 || isempty(rowRange)
        rows = 1:Ny;
    else
        rows = round(min(rowRange)):round(max(rowRange));
        rows = rows(rows >= 1 & rows <= Ny);
    end

    if isempty(rows)
        return
    end

    vals = residualMap(rows, :);
    vals = vals(isfinite(vals));
    if isempty(vals)
        return
    end

    medv = median(vals, 'omitnan');
end

% =========================================================================
% Convert real-world points to rotated-grid cell indices.
% =========================================================================
function [ix, iy] = local_points_to_cells(x, y, xi, yi, dx, dy, origin, e_s, e_n)
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
% Select one representative residual per A-point cell using the cell median.
% =========================================================================
function [selected_path_mask, cell_selected_idx, cell_median_val] = ...
    local_select_one_path_per_cell( ...
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

    for k = find(valid_mask)'
        xa = xA_path(k);
        ya = yA_path(k);
        vp = path_val(k);

        if ~isfinite(xa) || ~isfinite(ya) || ~isfinite(vp)
            continue
        end

        [ix_a, iy_a] = local_points_to_cells(xa, ya, xi, yi, dx, dy, origin, e_s, e_n);
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
% Local angle utilities mirroring the solver's residual-angle convention.
% =========================================================================
function [ang_out, shift] = local_wrapTo360_centerMedian(ang)
    shift = 0;
    ang_out = local_applyAngleShift(ang, shift);
end

function ang_out = local_applyAngleShift(ang, shift)
    if nargin < 2 || isempty(shift)
        shift = 0;
    end
    ang_out = mod(ang + shift + 180, 360) - 180;
end

% =========================================================================
% Helper for compact verbose printing of optional scalar values.
% =========================================================================
function x = local_printable_scalar(x)
    if isempty(x)
        x = NaN;
    end
end

% =========================================================================
% Find the row range that contains any active mask cell.
% =========================================================================
function [rowMin, rowMax] = local_active_row_range(evalMask)

    rowAny = any(evalMask, 2);
    rowIdx = find(rowAny);

    if isempty(rowIdx)
        rowMin = NaN;
        rowMax = NaN;
    else
        rowMin = min(rowIdx);
        rowMax = max(rowIdx);
    end
end