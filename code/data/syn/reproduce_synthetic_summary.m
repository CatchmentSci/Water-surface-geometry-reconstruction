function outputs = reproduce_synthetic_summary(archiveRoot, outputFolder)
%REPRODUCE_SYNTHETIC_SUMMARY Rebuild the accepted autocorrelation summary.
%   OUTPUTS = REPRODUCE_SYNTHETIC_SUMMARY(ARCHIVEROOT, OUTPUTFOLDER)
%   analyses the deposited synthetic solver checkpoints and writes the
%   autocorrelation summary used by Figure 5 and Appendix Table B1.
%
% The script scans files named like either of these equivalent formats:
%   syn_s1-5_solver_inputs_case_amp0.02_w0.63_ckpt.mat
%   syn_solver_inputs_0pt88_case_amp0.02_w0.63_ckpt.mat
%
% Both examples above are interpreted as the same synthetic input metadata:
%   h = 0.88 m, A = 0.02 m, lambda = 0.63 m, group S1-S5.
%
% Depth-coded files support h = 0.88, 1.10, 1.50, 1.70 and 1.79 m.
% Depth is parsed directly from the filename token, e.g.:
%   syn_solver_inputs_1pt50_case_amp0.02_w0.63_ckpt.mat -> h = 1.50 m
%   syn_solver_inputs_0pt88_case_amp0.08_w1.52_ckpt.mat -> h = 0.88 m
% Only depths with an explicit mapping are forced onto an S-group/case-number
% mapping. Other supported depth-coded files are analysed as unmapped cases.
%
% The prescribed amplitude, wavelength and depth are parsed from each filename.
% KLT_select_wave_checkpoint_refined supplies the reconstructed amplitude and
% median range-WSG autocorrelation wavelength and robust amplitude used by
% the published outputs. PMUSIC is deliberately disabled.
%
% The velocity difference uses the deep-water gravity-wave approximation:
%   U_s = sqrt(g*lambda/(2*pi))
%   Delta U [%] = 100*(U(lambda_est) - U(lambda_true))/U(lambda_true)
%
arguments
    archiveRoot (1, 1) string
    outputFolder (1, 1) string = fullfile(fileparts(mfilename('fullpath')), "output")
end

synFolder = fullfile(archiveRoot, "syn", "outputs");
if ~isfolder(outputFolder)
    mkdir(outputFolder);
end
outputAutocorrCsvFile = fullfile(outputFolder, ...
    "synthetic_wave_range_autocorr_summary.csv");
outputAutocorrMatFile = fullfile(outputFolder, ...
    "synthetic_wave_range_autocorr_summary.mat");

% Use the repository selector shared with the real-case data builders.
thisScriptDir = fileparts(mfilename('fullpath'));
selectorDir = fullfile(fileparts(thisScriptDir), "videos");
addpath(selectorDir, '-begin');
clear KLT_select_wave_checkpoint_refined
rehash;

fprintf('\nUsing selector: %s\n', which('KLT_select_wave_checkpoint_refined'));

% Selector options.
selectorOpts = struct;
selectorOpts.stabilityAssessmentMode = 'first_populated_10m_from_camera';
selectorOpts.firstSensedLength_m = 10;
selectorOpts.cameraRowOrder = 'ascending';
selectorOpts.displaySummary = false;

% The publication workflow uses autocorrelation only.
selectorOpts.usePmusicWavelength = false;

% Match the single-case caller's stop-map settings explicitly.
selectorOpts.requireForwardAmplitudePlateau = true;
selectorOpts.weakForwardWindow = 5;
selectorOpts.weakMaxForwardAmpGrowth = 0.05;
selectorOpts.weakMaxForwardAmpAbsGrowth_m = 0.003;
selectorOpts.strongForwardWindow = 5;
selectorOpts.strongMaxForwardAmpGrowth = 0.05;
selectorOpts.strongMaxForwardAmpAbsGrowth_m = 0.003;
selectorOpts.noPlateauFallbackMode = 'latest_supported';

%% ------------------------------------------------------------------------
% Find matching checkpoint files
%% ------------------------------------------------------------------------
if ~isfolder(synFolder)
    error('Syn folder does not exist: %s', synFolder);
end

% Support both historical group-coded filenames and newer depth-coded
% filenames. These two examples are treated as equivalent metadata:
%   syn_s1-5_solver_inputs_case_amp0.02_w0.63_ckpt.mat
%   syn_solver_inputs_0pt88_case_amp0.02_w0.63_ckpt.mat
patterns = { ...
    'syn_s*_solver_inputs_case_amp*_w*_ckpt.mat', ...
    'syn_solver_inputs_*_case_amp*_w*_ckpt.mat'};

allFiles = [];
for ip = 1:numel(patterns)
    thisFiles = dir(fullfile(synFolder, patterns{ip}));
    allFiles = [allFiles; thisFiles(:)]; %#ok<AGROW>
end
allFiles = unique_file_list(allFiles);

if isempty(allFiles)
    error(['No matching .mat files found in %s. Expected filenames like:\n', ...
           '  syn_s1-5_solver_inputs_case_amp0.02_w0.63_ckpt.mat\n', ...
           '  syn_solver_inputs_0pt88_case_amp0.02_w0.63_ckpt.mat\n', ...
           '  syn_solver_inputs_1pt50_case_amp0.02_w0.63_ckpt.mat'], synFolder);
end

% Keep only files that contain a checkpoint variable.
files = allFiles([]);
for ii = 1:numel(allFiles)
    fpath = fullfile(allFiles(ii).folder, allFiles(ii).name);
    vars = who('-file', fpath);
    if ismember('checkpoint', vars)
        files(end+1) = allFiles(ii); %#ok<SAGROW>
    end
end

if isempty(files)
    error('Matching files were found, but none contained a variable named checkpoint.');
end

%% ------------------------------------------------------------------------
% Run analysis over files
%% ------------------------------------------------------------------------
rows = repmat(empty_result_row(), numel(files), 1);
seenInputKeys = strings(0, 1);

for ii = 1:numel(files)
    fileName = files(ii).name;
    filePath = fullfile(files(ii).folder, fileName);

    meta = parse_synthetic_case_filename(fileName);
    if ~meta.valid
        warning('Skipping file because its name could not be parsed: %s', fileName);
        continue
    end

    if any(seenInputKeys == string(meta.inputKey))
        warning('Skipping duplicate synthetic input %s from file: %s', meta.inputKey, fileName);
        continue
    end
    seenInputKeys(end+1, 1) = string(meta.inputKey); %#ok<SAGROW>

    fprintf('\n====================================================================\n');
    fprintf('Analysing %s\n', fileName);
    fprintf('Parsed: case %s, h = %.3g m, A = %.5g m, lambda = %.5g m, format = %s\n', ...
        meta.caseLabelPlain, meta.depth_m, meta.amplitude_m, meta.wavelength_m, meta.fileNameFormat);

    % Populate the prescribed/input columns before running the selector.
    % If the selector fails, this row is still
    % retained in the output table and the calculated fields remain NaN.
    % The LaTeX writer renders those NaNs as '---'.
    lambdaTrue_m = meta.wavelength_m;
    Utrue_mps = deep_water_velocity_from_wavelength(lambdaTrue_m);

    rows(ii).caseNumber = meta.caseNumber;
    rows(ii).caseLabel = meta.caseLabelPlain;
    rows(ii).caseLabelLatex = meta.caseLabelLatex;
    rows(ii).fileName = string(fileName);
    rows(ii).groupLabel = string(meta.groupLabel);
    rows(ii).fileNameFormat = string(meta.fileNameFormat);
    rows(ii).inputKey = string(meta.inputKey);
    rows(ii).h_m = meta.depth_m;
    rows(ii).lambda_m = lambdaTrue_m;
    rows(ii).A_m = meta.amplitude_m;
    rows(ii).UtrueDeepWater_mps = Utrue_mps;
    rows(ii).analysisStatus = "pending";

    try
        S = load(filePath, 'checkpoint');

        [result, mapTable] = KLT_select_wave_checkpoint_refined(S.checkpoint, selectorOpts);

        [~, selectedMapRow] = resolve_stop_map_for_batch(result, mapTable);
        selectionReason = describe_selected_map_reason(result, mapTable, selectedMapRow);

        % Use the median range-WSG autocorrelation wavelength returned by
        % the selector; no second spectral calculation is required here.
        lambdaAutocorrEst_m = result.wavelength_m;
        UautocorrEst_mps = ...
            deep_water_velocity_from_wavelength(lambdaAutocorrEst_m);

        rows(ii).lambdaEstB_rangeWseAutocorr_m = lambdaAutocorrEst_m;
        rows(ii).velocityDiffB_rangeWseAutocorr_pct = ...
            percent_difference(UautocorrEst_mps, Utrue_mps);
        rows(ii).Aest_m = result.amplitude_m;
        rows(ii).selectedMapIndex = result.selectedMapIndex;
        rows(ii).selectedMapRow = selectedMapRow;
        rows(ii).stopMethod = string(result.stopMethod);
        rows(ii).selectionReason = string(selectionReason);
        rows(ii).nRangeWseAutocorrProfiles = result.nWavelengthProfiles;
        rows(ii).UestBAutocorrDeepWater_mps = UautocorrEst_mps;
        rows(ii).analysisStatus = "ok";

        suppressWaveAmpOutputs = should_suppress_wave_amplitude_outputs(result.stopMethod, selectionReason);
        rows(ii).suppressWaveAmpOutputs = suppressWaveAmpOutputs;
        if suppressWaveAmpOutputs
            rows(ii) = suppress_wave_amplitude_fields(rows(ii));
        end

    catch ME
        warning('Analysis failed for %s. Output row will be retained with calculated fields as --- in LaTeX. Reason: %s', ...
            fileName, ME.message);
        rows(ii).analysisStatus = "failed";
        rows(ii).failureMessage = string(ME.message);
        continue
    end

    fprintf('\n=== Batch summary for %s ===\n', meta.caseLabelPlain);
    if suppressWaveAmpOutputs
        fprintf('Estimated wavelength and amplitude: --- (latest_supported stop criterion)\n');
    else
        fprintf('Prescribed lambda: %.5g m\n', lambdaTrue_m);
        fprintf('Selector WSG autocorrelation wavelength median: %.5g m [IQR %.5g, %.5g], n = %d\n', ...
            lambdaAutocorrEst_m, result.wavelengthIQR_m(1), ...
            result.wavelengthIQR_m(2), result.nWavelengthProfiles);
        fprintf('Deep-water U true: %.5g m/s; autocorrelation U: %.5g m/s (%+.3f %%)\n', ...
            Utrue_mps, ...
            UautocorrEst_mps, ...
            rows(ii).velocityDiffB_rangeWseAutocorr_pct);
        fprintf('Prescribed A: %.5g m; estimated A: %.5g m\n', meta.amplitude_m, result.amplitude_m);
    end
    fprintf('Selected map reason: %s\n', selectionReason);
end

% Remove rows that failed to parse. Rows where analysis failed are retained.
validRows = arrayfun(@(r) strlength(r.fileName) > 0, rows);
rows = rows(validRows);

if isempty(rows)
    error('No valid synthetic checkpoint cases were analysed.');
end

T = struct2table(rows);
% Sort output table by hydraulic setup, not by inferred case number.
% Primary: depth h, then wavelength lambda, then amplitude A.
% caseNumber is retained only as a final tie-breaker for legacy S-labelled files.
T = sortrows(T, {'h_m','lambda_m','A_m','caseNumber'});

%% ------------------------------------------------------------------------
% Write outputs
%% ------------------------------------------------------------------------
writetable(T, outputAutocorrCsvFile);
save(outputAutocorrMatFile, 'T', 'selectorOpts', 'synFolder', ...
    'outputFolder');

fprintf('\n====================================================================\n');
fprintf('Batch analysis complete.\n');
fprintf('Autocorrelation CSV written to : %s\n', outputAutocorrCsvFile);
fprintf('Autocorrelation MAT written to : %s\n', outputAutocorrMatFile);
outputs = struct('csvFile', outputAutocorrCsvFile, ...
    'matFile', outputAutocorrMatFile, 'table', T, ...
    'checkpointFolder', synFolder);
end

%% ========================================================================
% Local helper functions
%% ========================================================================
function create_summary_relationship_figure( ...
        T, outputFigurePng, outputFigurePdf, ...
        relationshipFigureOpts, wavelengthMethod)
    if nargin < 4 || isempty(relationshipFigureOpts)
        relationshipFigureOpts = struct;
    end

    if nargin < 5 || isempty(wavelengthMethod)
        wavelengthMethod = 'pmusic';
    end
    methodSpec = wavelength_method_spec(wavelengthMethod);

    coolwarmFile = local_get_option( ...
        relationshipFigureOpts, 'coolwarmFile', '');

    if isempty(T) || height(T) == 0
        warning('Summary relationship figure was not created because T is empty.');
        return
    end

    if ismember('analysisStatus', T.Properties.VariableNames)
        okMask = string(T.analysisStatus) == "ok";
    else
        okMask = true(height(T), 1);
    end

    if ismember('suppressWaveAmpOutputs', T.Properties.VariableNames)
        okMask = okMask & ~logical(T.suppressWaveAmpOutputs);
    end

    lambdaTrue = T.lambda_m;
    lambdaEst = T.(methodSpec.wavelengthField);
    ampTrue = T.A_m;
    ampEst = T.Aest_m;

    validWave = okMask & isfinite(lambdaTrue) & isfinite(lambdaEst) & isfinite(ampTrue);
    validAmp = okMask & isfinite(ampTrue) & isfinite(ampEst) & isfinite(lambdaTrue);

    if ~any(validWave) && ~any(validAmp)
        warning('Summary relationship figure was not created because no finite analysed rows were available.');
        return
    end

    exportFigureWidthIn = 5.5;
    exportFigureHeightIn = 2.75;
    screenScale = local_get_option(relationshipFigureOpts, 'screenScale', 1);
    figureWidthIn = exportFigureWidthIn * screenScale;
    figureHeightIn = exportFigureHeightIn * screenScale;
    plotFontName = local_get_option(relationshipFigureOpts, 'fontName', 'Arial');
    axesFontSize = local_get_option(relationshipFigureOpts, 'axesFontSize', 8);
    labelFontSize = local_get_option(relationshipFigureOpts, 'labelFontSize', axesFontSize);
    colorbarFontSize = local_get_option(relationshipFigureOpts, 'colorbarFontSize', axesFontSize);
    colorbarLabelFontSize = local_get_option( ...
        relationshipFigureOpts, 'colorbarLabelFontSize', labelFontSize);
    panelLabelFontSize = local_get_option(relationshipFigureOpts, 'panelLabelFontSize', axesFontSize + 1);
    axesLineWidth = local_get_option(relationshipFigureOpts, 'axesLineWidth', 0.75);
    markerArea = local_get_option(relationshipFigureOpts, 'markerArea', 26);

    fig = figure('Color', 'w', 'Units', 'inches', ...
        'Position', [0.5 0.5 figureWidthIn figureHeightIn], ...
        'PaperUnits', 'inches', ...
        'PaperPosition', [0 0 exportFigureWidthIn exportFigureHeightIn], ...
        'PaperSize', [exportFigureWidthIn exportFigureHeightIn], ...
        'InvertHardcopy', 'off');

    tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    ax1 = nexttile(tl, 1);
    hold(ax1, 'on');
    grid(ax1, 'on');
    box(ax1, 'on');

    if any(validWave)
        scatter(ax1, lambdaTrue(validWave), lambdaEst(validWave), markerArea, ...
            ampTrue(validWave), 'o', 'filled', ...
            'MarkerFaceColor', 'flat', 'MarkerEdgeColor', [0.2 0.2 0.2], ...
            'LineWidth', 0.6);
    end

    colormap(ax1, local_external_colormap(coolwarmFile, 256));
    climA = local_color_limits(ampTrue(validWave));
    if ~isempty(climA)
        clim(ax1, climA);
    end
    cb1 = colorbar(ax1);
    cb1.Label.String = '$A$ (m)';
    cb1.Label.Interpreter = 'latex';
    cb1.FontName = plotFontName;
    cb1.FontSize = colorbarFontSize;
    cb1.Label.FontName = plotFontName;
    cb1.Label.FontSize = colorbarLabelFontSize;
    if ~isempty(climA)
        % Fix both tick positions and labels so tiled-layout resizing cannot
        % leave stale labels attached to newly calculated tick positions.
        cb1.Ticks = linspace(climA(1), climA(2), 5);
        cb1.TickLabels = compose('%.2f', cb1.Ticks);
    end

    waveVals = [lambdaTrue(validWave); lambdaEst(validWave)];
    limsWave = make_equal_limits(waveVals);
    if ~isempty(limsWave)
        plot(ax1, limsWave, limsWave, 'k--', 'LineWidth', 0.9);
        xlim(ax1, limsWave);
        ylim(ax1, limsWave);
        waveTicks = ceil(limsWave(1)):floor(limsWave(2));
        if ~isempty(waveTicks)
            xticks(ax1, waveTicks);
            yticks(ax1, waveTicks);
        end
    end

    xlabel(ax1, '$\lambda$ (m)', 'Interpreter', 'latex');
    ylabel(ax1, '$\lambda_{\mathrm{est}}$ (m)', 'Interpreter', 'latex');
    ax1.XLabel.FontName = plotFontName;
    ax1.XLabel.FontSize = labelFontSize;
    ax1.YLabel.FontName = plotFontName;
    ax1.YLabel.FontSize = labelFontSize;
    axis(ax1, 'square');
    add_panel_label_local(ax1, '(a)');

    ax2 = nexttile(tl, 2);
    hold(ax2, 'on');
    grid(ax2, 'on');
    box(ax2, 'on');

    if any(validAmp)
        scatter(ax2, ampTrue(validAmp), ampEst(validAmp), markerArea, ...
            lambdaTrue(validAmp), 'o', 'filled', ...
            'MarkerFaceColor', 'flat', 'MarkerEdgeColor', [0.2 0.2 0.2], ...
            'LineWidth', 0.6);
    end
    colormap(ax2, local_external_colormap(coolwarmFile, 256));
    climB = local_color_limits(lambdaTrue(validAmp));
    if ~isempty(climB)
        clim(ax2, climB);
    end
    cb2 = colorbar(ax2);
    cb2.Label.String = '$\lambda$ (m)';
    cb2.Label.Interpreter = 'latex';
    cb2.FontName = plotFontName;
    cb2.FontSize = colorbarFontSize;
    cb2.Label.FontName = plotFontName;
    cb2.Label.FontSize = colorbarLabelFontSize;

    ampVals = [ampTrue(validAmp); ampEst(validAmp)];
    limsAmp = make_equal_limits(ampVals);
    if ~isempty(limsAmp)
        plot(ax2, limsAmp, limsAmp, 'k--', 'LineWidth', 0.9);
        xlim(ax2, limsAmp);
        ylim(ax2, limsAmp);
        ampTickStep = 0.02;
        ampTickIndices = ceil(limsAmp(1) / ampTickStep):floor(limsAmp(2) / ampTickStep);
        ampTicks = ampTickIndices * ampTickStep;
        if ~isempty(ampTicks)
            xticks(ax2, ampTicks);
            yticks(ax2, ampTicks);
            xtickformat(ax2, '%.2f');
            ytickformat(ax2, '%.2f');
        end
    end

    xlabel(ax2, '$A$ (m)', 'Interpreter', 'latex');
    ylabel(ax2, '$A_{\mathrm{est}}$ (m)', 'Interpreter', 'latex');
    ax2.XLabel.FontName = plotFontName;
    ax2.XLabel.FontSize = labelFontSize;
    ax2.YLabel.FontName = plotFontName;
    ax2.YLabel.FontSize = labelFontSize;
    axis(ax2, 'square');
    add_panel_label_local(ax2, '(b)');

    set([ax1 ax2], 'FontName', plotFontName, 'FontSize', axesFontSize, ...
        'LineWidth', axesLineWidth, 'TickDir', 'out', 'Layer', 'top', ...
        'GridAlpha', 0.12);
    panelLabels = findall(fig, 'Tag', 'panelLabel');
    set(panelLabels, 'FontName', plotFontName, 'FontSize', panelLabelFontSize);

    drawnow;
    exportgraphics(fig, outputFigurePng, 'Resolution', 600);
    exportgraphics(fig, outputFigurePdf, 'ContentType', 'vector');
end

function value = local_get_option(S, fieldName, defaultValue)
    if isstruct(S) && isfield(S, fieldName) && ~isempty(S.(fieldName))
        value = S.(fieldName);
    else
        value = defaultValue;
    end
end

function climVals = local_color_limits(vals)
    vals = vals(isfinite(vals));
    if isempty(vals)
        climVals = [];
        return
    end
    vmin = min(vals);
    vmax = max(vals);
    if vmin == vmax
        pad = max(abs(vmin) * 0.05, 1e-6);
        climVals = [vmin - pad, vmax + pad];
    else
        climVals = [vmin, vmax];
    end
end

function cmap = local_external_colormap(cmapFile, n)
    if nargin < 2 || isempty(n)
        n = 256;
    end

    cmap = parula(n);

    if isempty(cmapFile) || ~isfile(cmapFile)
        warning('Could not find external colormap file: %s. Using parula instead.', cmapFile);
        return
    end

    try
        txt = fileread(cmapFile);

        % Parse the numeric block assigned to c = [ ... ];
        % (?s) enables dot to match newlines.
        token = regexp(txt, '(?s)c\s*=\s*\[(.*?)\];', 'tokens', 'once');
        if isempty(token)
            warning('Could not parse colormap data from %s. Using parula instead.', cmapFile);
            return
        end

        rawBlock = token{1};

        % Split into rows using semicolons or line breaks
        rowTokens = regexp(rawBlock, '[;\r\n]+', 'split');
        rowTokens = rowTokens(~cellfun(@isempty, strtrim(rowTokens)));

        c = zeros(0,3);
        for iRow = 1:numel(rowTokens)
            vals = sscanf(rowTokens{iRow}, '%f').';
            if isempty(vals)
                continue
            end
            if numel(vals) ~= 3
                warning('Row %d in %s does not contain exactly 3 values. Using parula instead.', ...
                    iRow, cmapFile);
                cmap = parula(n);
                return
            end
            c(end+1, :) = vals; %#ok<AGROW>
        end

        if isempty(c) || size(c,2) ~= 3
            warning('Parsed colormap data from %s is invalid. Using parula instead.', cmapFile);
            cmap = parula(n);
            return
        end

        % Normalize if stored as 0-255 RGB
        if max(c(:)) > 1
            c = c / 255;
        end

        nBase = size(c,1);

        if nBase == 1
            cmap = repmat(c, n, 1);
            return
        end

        if n == 1
            cmap = c(round((nBase + 1)/2), :);
            return
        end

        xBase  = linspace(1, nBase, nBase);
        xQuery = linspace(1, nBase, n);

        r = interp1(xBase, c(:,1), xQuery, 'linear');
        g = interp1(xBase, c(:,2), xQuery, 'linear');
        b = interp1(xBase, c(:,3), xQuery, 'linear');

        cmap = [r(:) g(:) b(:)];
        cmap = max(0, min(1, cmap));

    catch ME
        warning('Failed to read external colormap file %s (%s). Using parula instead.', cmapFile, ME.message);
        cmap = parula(n);
    end
end

function lims = make_equal_limits(vals)
    vals = vals(isfinite(vals));
    if isempty(vals)
        lims = [];
        return
    end

    vmin = min(vals);
    vmax = max(vals);
    if ~isfinite(vmin) || ~isfinite(vmax)
        lims = [];
        return
    end

    if vmin == vmax
        pad = max(abs(vmin) * 0.075, 1e-3);
        lims = [vmin - pad, vmax + pad];
    else
        pad = 0.075 * (vmax - vmin);
        lims = [vmin - pad, vmax + pad];
    end
end

function add_panel_label_local(ax, labelText)
    text(ax, 0.02, 0.97, labelText, 'Units', 'normalized', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
        'FontWeight', 'bold', 'Tag', 'panelLabel');
end

function R = empty_result_row()
    R = struct;
    R.caseNumber = NaN;
    R.caseLabel = "";
    R.caseLabelLatex = "";
    R.fileName = "";
    R.groupLabel = "";
    R.fileNameFormat = "";
    R.inputKey = "";
    R.h_m = NaN;
    R.lambda_m = NaN;
    R.lambdaEstB_rangeWsePmusic_m = NaN;
    R.velocityDiffB_rangeWsePmusic_pct = NaN;
    R.lambdaEstB_rangeWseAutocorr_m = NaN;
    R.velocityDiffB_rangeWseAutocorr_pct = NaN;
    R.A_m = NaN;
    R.Aest_m = NaN;
    R.selectedMapIndex = NaN;
    R.selectedMapRow = NaN;
    R.stopMethod = "";
    R.selectionReason = "";
    R.analysisStatus = "";
    R.failureMessage = "";
    R.nRangeWsePmusicProfiles = NaN;
    R.nRangeWseAutocorrProfiles = NaN;
    R.UtrueDeepWater_mps = NaN;
    R.UestBDeepWater_mps = NaN;
    R.UestBAutocorrDeepWater_mps = NaN;
    R.suppressWaveAmpOutputs = false;
end

function meta = parse_synthetic_case_filename(fileName)
    [~, nameOnly, ~] = fileparts(fileName);

    % Supported filename formats:
    %   1) syn_s1-5_solver_inputs_case_amp0.02_w0.63_ckpt.mat
    %   2) syn_solver_inputs_0pt88_case_amp0.02_w0.63_ckpt.mat
    %
    % In format 1, depth is inferred from the explicit S-group.
    % In format 2, depth is read directly from the depth token. Supported
    % depth-coded filenames currently include:
    %   0.88, 1.10, 1.50, 1.70 and 1.79 m
    % Only depths with an explicit mapping are converted to S-groups / case
    % numbers. Other supported depths are analysed and retained in the
    % output table, but left as unmapped depth-coded cases.
    numberExpr = '\d+(?:(?:\.|pt)\d+)?';

    exprGroup = ['^syn_s(?<sStart>\d+)-(?<sEnd>\d+)_solver_inputs_case_', ...
                 'amp(?<amp>', numberExpr, ')_w(?<w>', numberExpr, ')_ckpt$'];
    exprDepth = ['^syn_solver_inputs_(?<depth>', numberExpr, ')_case_', ...
                 'amp(?<amp>', numberExpr, ')_w(?<w>', numberExpr, ')_ckpt$'];

    tokGroup = regexp(nameOnly, exprGroup, 'names', 'once');
    tokDepth = regexp(nameOnly, exprDepth, 'names', 'once');

    meta = struct;
    meta.valid = false;
    meta.sStart = NaN;
    meta.sEnd = NaN;
    meta.groupLabel = '';
    meta.caseNumber = NaN;
    meta.caseLabelPlain = '';
    meta.caseLabelLatex = '';
    meta.depth_m = NaN;
    meta.amplitude_m = NaN;
    meta.wavelength_m = NaN;
    meta.fileNameFormat = '';
    meta.inputKey = '';

    if ~isempty(tokGroup)
        meta.fileNameFormat = 'group-coded';
        meta.sStart = str2double(tokGroup.sStart);
        meta.sEnd = str2double(tokGroup.sEnd);
        meta.depth_m = depth_from_group(meta.sStart, meta.sEnd);
        meta.amplitude_m = parse_filename_number(tokGroup.amp);
        meta.wavelength_m = parse_filename_number(tokGroup.w);

    elseif ~isempty(tokDepth)
        meta.fileNameFormat = 'depth-coded';
        meta.depth_m = parse_filename_number(tokDepth.depth);
        meta.amplitude_m = parse_filename_number(tokDepth.amp);
        meta.wavelength_m = parse_filename_number(tokDepth.w);

        if ~is_supported_depth_token(meta.depth_m)
            warning('Skipping depth-coded file with unsupported depth %.5g m: %s', meta.depth_m, fileName);
            return
        end

        % Only map depth-coded files to S-groups where a mapping is defined.
        % h = 1.10, 1.70 and 1.79 m are currently analysed but intentionally
        % left unmapped.
        [meta.sStart, meta.sEnd] = group_from_depth(meta.depth_m);

    else
        return
    end

    % A valid row only requires depth, amplitude and wavelength. S-group and
    % case number are optional because some supported depth-coded cases
    % (e.g. 1.10, 1.70, 1.79 m) have no current mapping.
    meta.valid = isfinite(meta.depth_m) && isfinite(meta.amplitude_m) && ...
        isfinite(meta.wavelength_m);

    if ~meta.valid
        return
    end

    meta.inputKey = sprintf('h%.5g_A%.5g_w%.5g', meta.depth_m, meta.amplitude_m, meta.wavelength_m);

    if isfinite(meta.sStart) && isfinite(meta.sEnd)
        meta.groupLabel = sprintf('syn_s%d-%d', meta.sStart, meta.sEnd);
        meta.caseNumber = infer_case_number(meta.sStart, meta.sEnd, meta.amplitude_m, meta.wavelength_m);
    else
        meta.groupLabel = sprintf('h%.2f_unmapped', meta.depth_m);
        meta.caseNumber = NaN;
    end

    if isfinite(meta.caseNumber)
        meta.caseLabelPlain = sprintf('S%d', meta.caseNumber);
        if meta.caseNumber < 10
            meta.caseLabelLatex = sprintf('S$_%d$', meta.caseNumber);
        else
            meta.caseLabelLatex = sprintf('S$_{%d}$', meta.caseNumber);
        end
    else
        % No S-case mapping exists for this input, so keep the label explicit
        % rather than inventing an S number.
        meta.caseLabelPlain = sprintf('h%.2f_A%.5g_w%.5g', meta.depth_m, meta.amplitude_m, meta.wavelength_m);
        meta.caseLabelLatex = sprintf('$h=%.2f$', meta.depth_m);
    end
end

function h = depth_from_group(sStart, sEnd)
    if sStart == 1 && sEnd == 5
        h = 0.88;
    elseif sStart == 6 && sEnd == 10
        h = 1.50;
    elseif sStart == 11 && sEnd == 15
        % The old 2.16 m depth is no longer used. Group-coded S11-S15 files,
        % if present, are treated as the current high-depth set.
        h = 1.79;
    else
        h = NaN;
    end
end

function tf = is_supported_depth_token(depth_m)
    tol = 5e-3;
    tf = any(abs(depth_m - [0.88 1.10 1.50 1.70 1.79]) <= tol);
end

function [sStart, sEnd] = group_from_depth(depth_m)
    sStart = NaN;
    sEnd = NaN;

    tol = 5e-3;
    if abs(depth_m - 0.88) <= tol
        sStart = 1;
        sEnd = 5;
    elseif abs(depth_m - 1.50) <= tol
        sStart = 6;
        sEnd = 10;
    elseif abs(depth_m - 1.10) <= tol
        % Supported depth-coded case, but no S-group/case-number mapping defined.
        sStart = NaN;
        sEnd = NaN;
    elseif abs(depth_m - 1.70) <= tol
        % Supported depth-coded case, but no S-group/case-number mapping defined.
        sStart = NaN;
        sEnd = NaN;
    elseif abs(depth_m - 1.79) <= tol
        % Supported depth-coded case, but no S-group/case-number mapping defined.
        sStart = NaN;
        sEnd = NaN;
    end
end

function x = parse_filename_number(token)
    token = char(token);
    token = strrep(token, 'pt', '.');
    x = str2double(token);
end

function files = unique_file_list(files)
    if isempty(files)
        return
    end

    fullNames = strings(numel(files), 1);
    for ii = 1:numel(files)
        fullNames(ii) = string(fullfile(files(ii).folder, files(ii).name));
    end

    [~, keepIdx] = unique(fullNames, 'stable');
    files = files(keepIdx);
end

function caseNum = infer_case_number(sStart, sEnd, amp, w)
    caseNum = NaN;
    tol = 5e-3;

    if sStart == 1 && sEnd == 5
        known = [
            1, 0.02, 0.63
            2, 0.02, 1.52
            3, 0.08, 1.52
            4, 0.02, 3.46
            5, 0.10, 3.46];
    elseif sStart == 6 && sEnd == 10
        known = [
            6, 0.02, 0.63
            7, 0.02, 1.58
            8, 0.08, 1.58
            9, 0.02, 3.51
            10, 0.10, 3.51];
    elseif sStart == 11 && sEnd == 15
        known = [
            11, 0.02, 0.63
            12, 0.02, 1.62
            13, 0.08, 1.62
            14, 0.02, 3.62
            15, 0.10, 3.62];
    else
        known = [];
    end

    if isempty(known)
        return
    end

    d = abs(known(:,2) - amp) + abs(known(:,3) - w);
    [dmin, idx] = min(d);
    if dmin <= 2*tol
        caseNum = known(idx, 1);
    end
end

function selectionReason = describe_selected_map_reason(result, mapTable, selectedMapRow)
    selectionReason = "";

    if ~isempty(mapTable) && istable(mapTable) && isfinite(selectedMapRow) && ...
            selectedMapRow >= 1 && selectedMapRow <= height(mapTable)

        preferredCols = {'stopDiagnostics','selectionReason','stopReason', ...
                         'reason','selectionMethod','stopMethod'};

        for kk = 1:numel(preferredCols)
            vn = preferredCols{kk};
            if ismember(vn, mapTable.Properties.VariableNames)
                val = mapTable{selectedMapRow, vn};

                if iscell(val)
                    val = val{1};
                end

                if isstring(val)
                    val = val(1);
                end

                if ischar(val)
                    val = string(val);
                elseif isnumeric(val) && isscalar(val) && isfinite(val)
                    val = string(num2str(val));
                end

                if isstring(val) && strlength(strtrim(val)) > 0
                    selectionReason = strtrim(val);
                    break
                end
            end
        end
    end

    if strlength(selectionReason) == 0
        if isfield(result, 'stopMethod') && strlength(string(result.stopMethod)) > 0
            txt = string(result.stopMethod);
            txt = replace(txt, "_", " ");
            selectionReason = txt;
        else
            selectionReason = "selected map reason unavailable";
        end
    end

    mapTxt = "";
    rowTxt = "";

    if isfield(result, 'selectedMapIndex') && isfinite(result.selectedMapIndex)
        mapTxt = "map " + string(round(result.selectedMapIndex));
    end

    if isfinite(selectedMapRow)
        rowTxt = "row " + string(round(selectedMapRow));
    end

    if strlength(mapTxt) > 0 && strlength(rowTxt) > 0
        selectionReason = selectionReason + " (" + mapTxt + ", " + rowTxt + ")";
    elseif strlength(mapTxt) > 0
        selectionReason = selectionReason + " (" + mapTxt + ")";
    elseif strlength(rowTxt) > 0
        selectionReason = selectionReason + " (" + rowTxt + ")";
    end
end


function [stopMap, selectedMapRow] = resolve_stop_map_for_batch(result, mapTable)
    selectedMapRow = NaN;
    if isfield(result, 'selectedMapRow') && isfinite(result.selectedMapRow)
        selectedMapRow = round(result.selectedMapRow);
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
            selectedMapRow = exactRow;
            stopMap = selectedValue;
            return
        end
        if selectedValue >= 1 && selectedValue <= height(mapTable)
            selectedMapRow = selectedValue;
            stopMap = round(mapTable.mapIndex(selectedMapRow));
            return
        end
    end
end


function U = deep_water_velocity_from_wavelength(lambda_m)
    g = 9.81;
    if ~isfinite(lambda_m) || lambda_m <= 0
        U = NaN;
    else
        U = sqrt(g .* lambda_m ./ (2*pi));
    end
end

function pct = percent_difference(estimatedValue, referenceValue)
    if ~isfinite(estimatedValue) || ~isfinite(referenceValue) || referenceValue == 0
        pct = NaN;
    else
        pct = 100 .* (estimatedValue - referenceValue) ./ referenceValue;
    end
end

function write_latex_wave_summary_table(T, outTexFile, wavelengthMethod)
    if nargin < 3 || isempty(wavelengthMethod)
        wavelengthMethod = 'pmusic';
    end
    methodSpec = wavelength_method_spec(wavelengthMethod);
    lambdaEstimate = T.(methodSpec.wavelengthField);
    velocityDifference = T.(methodSpec.velocityDifferenceField);

    fid = fopen(outTexFile, 'w');
    if fid < 0
        error('Could not open LaTeX output file for writing: %s', outTexFile);
    end
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, ['%% Auto-generated by ', ...
        'calling_syn_batch_analysis_latex_table_simplified.m\n']);
    fprintf(fid, '%% Requires \\usepackage{booktabs}\n');
    fprintf(fid, '\\begin{table}[htbp]\n');
    fprintf(fid, '\\centering\n');
    fprintf(fid, '\\setlength{\\tabcolsep}{8pt}\n');

    fprintf(fid, '\\begin{tabular}{c c c c c c c}\n');
    fprintf(fid, '\\toprule\n');
    fprintf(fid, 'Case & $h$~[\\mathrm{m}] & $\\lambda$~[\\mathrm{m}] & $\\lambda_{\\mathrm{est}}$~[\\mathrm{m}] & $\\Delta U$~[\\%%] & $A$~[\\mathrm{m}] & $A_{\\mathrm{est}}$~[\\mathrm{m}] \\\\\n');
    fprintf(fid, '\\midrule\n');

    for ii = 1:height(T)
        caseLabelLatexForTable = sequential_latex_case_label(ii);

        fprintf(fid, '%s & %s & %s & %s & %s & %s & %s \\\\\n', ...
            caseLabelLatexForTable, ...
            latex_num(T.h_m(ii), '%.2f'), ...
            latex_num(T.lambda_m(ii), '%.3f'), ...
            latex_num(lambdaEstimate(ii), '%.3f'), ...
            latex_num(velocityDifference(ii), '%+.2f'), ...
            latex_num(T.A_m(ii), '%.3f'), ...
            latex_num(T.Aest_m(ii), '%.3f'));
    end

    fprintf(fid, '\\bottomrule\n');
    fprintf(fid, '\\end{tabular}\n');
    fprintf(fid, '\\bigskip\n');

    fprintf(fid, ['\\caption{Characteristics associated with each numerical simulation, where water depth $h$, ', ...
                  'wavelength $\\lambda$, and amplitude $A$ are varied. ', ...
                  '$\\lambda_{\\mathrm{est}}$ is the median range WSE %s wavelength estimate. ', ...
                  '$\\Delta U$ is the signed percentage difference between the deep-water gravity-wave velocity ', ...
                  'calculated from the estimated wavelength and that calculated from the prescribed wavelength ', ...
                  'using $U_s=\\sqrt{g\\lambda/(2\\pi)}$. Estimated values are calculated as medians of streamwise ', ...
                  'transects with cross-stream spacing of 0.5~m.}\n'], ...
                  methodSpec.captionName);

    fprintf(fid, '\\label{%s}\n', methodSpec.tableLabel);
    fprintf(fid, '\\end{table}\n');
end

function label = sequential_latex_case_label(rowNumber)
    rowNumber = round(rowNumber);
    if rowNumber < 10
        label = sprintf('S$_%d$', rowNumber);
    else
        label = sprintf('S$_{%d}$', rowNumber);
    end
end

function s = latex_num(x, fmt)
    if isfinite(x)
        s = sprintf(fmt, x);
    else
        s = '---';
    end
end

function tf = should_suppress_wave_amplitude_outputs(stopMethod, selectionReason)
    stopMethod = lower(string(stopMethod));
    selectionReason = lower(string(selectionReason));

    tf = any(contains(stopMethod, "latest_supported")) || ...
         any(contains(selectionReason, "latest quality-supported")) || ...
         any(contains(selectionReason, "latest finite map"));
end

function R = suppress_wave_amplitude_fields(R)
    R.lambdaEstB_rangeWsePmusic_m = NaN;
    R.velocityDiffB_rangeWsePmusic_pct = NaN;
    R.lambdaEstB_rangeWseAutocorr_m = NaN;
    R.velocityDiffB_rangeWseAutocorr_pct = NaN;
    R.Aest_m = NaN;
    R.UestBDeepWater_mps = NaN;
    R.UestBAutocorrDeepWater_mps = NaN;
end

function spec = wavelength_method_spec(wavelengthMethod)
    switch lower(string(wavelengthMethod))
        case "pmusic"
            spec.wavelengthField = 'lambdaEstB_rangeWsePmusic_m';
            spec.velocityDifferenceField = ...
                'velocityDiffB_rangeWsePmusic_pct';
            spec.captionName = 'PMUSIC';
            spec.tableLabel = 'Table:syn_res';
        case {"autocorr", "autocorrelation"}
            spec.wavelengthField = 'lambdaEstB_rangeWseAutocorr_m';
            spec.velocityDifferenceField = ...
                'velocityDiffB_rangeWseAutocorr_pct';
            spec.captionName = 'autocorrelation';
            spec.tableLabel = 'Table:syn_res_autocorr';
        otherwise
            error('Unsupported wavelength method: %s', ...
                char(string(wavelengthMethod)));
    end
end
