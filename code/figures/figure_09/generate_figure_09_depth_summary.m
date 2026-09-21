function outputs = generate_figure_09_depth_summary(sourceFile, outputFolder)
%GENERATE_FIGURE_09_DEPTH_SUMMARY Rebuild the initial-velocity depth data.
%   Every admissible subcritical root in the 0.05--5.00 m interval is
%   identified by bracketing sign changes on a 0.01 m search grid. The
%   largest root is selected when multiple roots exist; the sole root is
%   retained when only one admissible root exists.

arguments
    sourceFile (1,1) string = fullfile(fileparts(mfilename("fullpath")), ...
        "..", "figure_08", "data", ...
        "wse_autocorrelation_velocity_method_sensitivity_summary.mat")
    outputFolder (1,1) string = fullfile(fileparts(mfilename("fullpath")), "data")
end

sourceFile = string(sourceFile);
outputFolder = string(outputFolder);
if ~isfile(sourceFile)
    error("Figure09Depth:MissingInput", ...
        "Required Figure 8 transect data were not found:\n%s", sourceFile);
end
if ~isfolder(outputFolder)
    mkdir(outputFolder);
end

source = load(sourceFile, "autocorrCaseTransectTables", "wavelengthFilter");
if ~isfield(source, "autocorrCaseTransectTables") || ...
        ~iscell(source.autocorrCaseTransectTables) || ...
        numel(source.autocorrCaseTransectTables) ~= 13
    error("Figure09Depth:InvalidInput", ...
        "Expected 13 ordered R1--R13 transect tables.");
end
validate_wavelength_filter(source);

cfg = inversion_config();
methods = ["uniform", "linear", "power"];
caseRows = repmat(empty_case_row(methods), 13, 1);
auditRows = cell(13 .* numel(methods), 1);
auditIndex = 0;

for ii = 1:13
    T = source.autocorrCaseTransectTables{ii};
    required = ["caseLabel", "discharge_Q_m3ps", "profileIndex", ...
        "crossSectionDepth_m", "acceptedAutocorrWavelength_m", ...
        "initialVelocityTrackedMedian_mps"];
    missing = setdiff(required, string(T.Properties.VariableNames));
    if ~isempty(missing)
        error("Figure09Depth:InvalidInput", ...
            "Case %d is missing field(s): %s", ii, strjoin(missing, ", "));
    end

    expectedLabel = "R" + ii;
    labels = unique(string(T.caseLabel));
    if numel(labels) ~= 1 || labels ~= expectedLabel
        error("Figure09Depth:UnexpectedCases", ...
            "Case table %d is not labelled %s.", ii, expectedLabel);
    end

    wavelength = double(T.acceptedAutocorrWavelength_m);
    velocity = double(T.initialVelocityTrackedMedian_mps);
    observedDepth = double(T.crossSectionDepth_m);
    validBase = isfinite(wavelength) & wavelength > 0 & ...
        isfinite(velocity) & velocity > 0 & ...
        isfinite(observedDepth) & observedDepth > 0;

    row = empty_case_row(methods);
    row.caseNumber = ii;
    row.caseLabel = expectedLabel;
    row.discharge_Q_m3ps = median(double(T.discharge_Q_m3ps), "omitnan");
    row.includedInFigure = ii <= 11;
    if ~row.includedInFigure
        row.exclusionReason = ...
            "Rejected erroneous autocorrelation wavelength reconstruction";
    end
    row.nTransects = height(T);
    row.nSourceAcceptedTransects = nnz(validBase);

    for mm = 1:numel(methods)
        method = methods(mm);
        shallow = nan(height(T), 1);
        selected = nan(height(T), 1);
        rootCount = zeros(height(T), 1);
        allRoots = cell(height(T), 1);

        for jj = 1:height(T)
            if ~validBase(jj)
                continue
            end
            roots = find_all_depth_roots( ...
                2*pi./wavelength(jj), velocity(jj), method, cfg);
            rootCount(jj) = numel(roots);
            allRoots{jj} = roots;
            if ~isempty(roots)
                shallow(jj) = roots(1);
                selected(jj) = roots(end);
            end
        end

        row.(method + "_nOneSolution") = nnz(validBase & rootCount == 1);
        row.(method + "_nTwoSolutions") = nnz(validBase & rootCount >= 2);
        row = add_pair_summary(row, observedDepth, shallow, ...
            validBase & isfinite(shallow), method + "_shallow");
        row = add_pair_summary(row, observedDepth, selected, ...
            validBase & isfinite(selected), method + "_deep");

        auditIndex = auditIndex + 1;
        auditRows{auditIndex} = table( ...
            repmat(ii, height(T), 1), repmat(expectedLabel, height(T), 1), ...
            double(T.profileIndex), repmat(method, height(T), 1), ...
            observedDepth, wavelength, velocity, rootCount, shallow, selected, ...
            strings(height(T), 1), ...
            VariableNames=["caseNumber", "caseLabel", "profileIndex", ...
            "velocityProfile", "observedDepth_m", "wavelength_m", ...
            "initialPlanarVelocity_mps", "nAdmissibleRoots", ...
            "shallowOrSoleRoot_m", "selectedDeepOrSoleRoot_m", "allRoots_m"]);
        for jj = 1:height(T)
            auditRows{auditIndex}.allRoots_m(jj) = ...
                strjoin(compose("%.15g", allRoots{jj}), ";");
        end
    end
    caseRows(ii) = row;
end

summaryTable = struct2table(caseRows);
auditTable = vertcat(auditRows{1:auditIndex});
summaryFile = fullfile(outputFolder, ...
    "wse_autocorrelation_uniform_linear_power_depth_summary.csv");
auditFile = fullfile(outputFolder, ...
    "wse_autocorrelation_depth_root_audit.csv");
writetable(summaryTable, summaryFile);
writetable(auditTable, auditFile);

outputs = struct("sourceFile", sourceFile, "summaryFile", summaryFile, ...
    "auditFile", auditFile, "summaryTable", summaryTable, ...
    "auditTable", auditTable);
end

function cfg = inversion_config()
cfg.g_mps2 = 9.81;
cfg.surfaceTension_Npm = 72.75e-3;
cfg.density_kgpm3 = 998.2;
cfg.alpha = 0.85;
cfg.depthSearchGrid_m = (0.01:0.01:10).';
cfg.reasonableDepthRange_m = [0.05 5.00];
cfg.depthRootResidualTolerance = 1e-8;
cfg.tangentRootResidualTolerance = 1e-4;
cfg.depthRootMergeTolerance_m = 1e-4;
end

function roots = find_all_depth_roots(k, velocity, method, cfg)
grid = cfg.depthSearchGrid_m;
residualFunction = @(h) profile_velocity_squared(k, h, method, cfg) - velocity.^2;
residual = arrayfun(residualFunction, grid);
finite = isfinite(residual);
roots = grid(finite & abs(residual) <= cfg.depthRootResidualTolerance);

crossings = find(finite(1:end-1) & finite(2:end) & ...
    residual(1:end-1).*residual(2:end) < 0);
for jj = crossings.'
    try
        roots(end+1, 1) = fzero(residualFunction, ...
            [grid(jj), grid(jj+1)]); %#ok<AGROW>
    catch
    end
end

for jj = 2:numel(grid)-1
    if finite(jj-1) && finite(jj) && finite(jj+1) && ...
            abs(residual(jj)) <= abs(residual(jj-1)) && ...
            abs(residual(jj)) <= abs(residual(jj+1))
        try
            candidate = fminbnd(@(h) abs(residualFunction(h)), ...
                grid(jj-1), grid(jj+1));
            if abs(residualFunction(candidate)) <= ...
                    cfg.tangentRootResidualTolerance
                roots(end+1, 1) = candidate; %#ok<AGROW>
            end
        catch
        end
    end
end

roots = roots(isfinite(roots) & roots >= cfg.reasonableDepthRange_m(1) & ...
    roots <= cfg.reasonableDepthRange_m(2));
roots = roots(velocity.^2 ./ (cfg.g_mps2 .* roots) <= 1 + 1e-10);
roots = sort(roots);
if ~isempty(roots)
    roots = roots([true; diff(roots) > cfg.depthRootMergeTolerance_m]);
end
end

function value = profile_velocity_squared(k, depth, method, cfg)
kh = k .* depth;
capillaryFactor = 1 + cfg.surfaceTension_Npm .* k.^2 ./ ...
    (cfg.density_kgpm3 .* cfg.g_mps2);
switch method
    case "uniform"
        value = cfg.g_mps2 .* depth .* capillaryFactor .* tanh(kh) ./ kh;
    case "linear"
        m = 2 .* (1-cfg.alpha);
        value = cfg.g_mps2 .* depth .* capillaryFactor .* tanh(kh) ./ ...
            (kh-m.*tanh(kh));
    case "power"
        n = 1./cfg.alpha - 1;
        s = sign(0.5-n);
        if s == 0, s = 1; end
        ratio = besseli(s.*(0.5-n), kh, 1) ./ ...
            besseli(-s.*(0.5+n), kh, 1);
        value = cfg.g_mps2 .* depth .* capillaryFactor .* ratio ./ kh;
    otherwise
        error("Figure09Depth:UnknownMethod", "Unknown method: %s", method);
end
end

function row = empty_case_row(methods)
row = struct("caseNumber", NaN, "caseLabel", "", ...
    "discharge_Q_m3ps", NaN, "includedInFigure", false, ...
    "exclusionReason", "", "nTransects", 0, ...
    "nSourceAcceptedTransects", 0);
for method = methods
    row.(method + "_nOneSolution") = 0;
    row.(method + "_nTwoSolutions") = 0;
    for branch = ["shallow", "deep"]
        prefix = method + "_" + branch;
        row.(prefix + "_nPairs") = 0;
        for suffix = ["_xMedian", "_xQ25", "_xQ75", "_xErrLow", ...
                "_xErrHigh", "_yMedian", "_yQ25", "_yQ75", ...
                "_yErrLow", "_yErrHigh"]
            row.(prefix + suffix) = NaN;
        end
    end
end
end

function row = add_pair_summary(row, x, y, valid, prefix)
x = double(x(valid));
y = double(y(valid));
row.(prefix + "_nPairs") = numel(x);
if isempty(x)
    return
end
xMedian = median(x); xQ25 = prctile(x, 25); xQ75 = prctile(x, 75);
yMedian = median(y); yQ25 = prctile(y, 25); yQ75 = prctile(y, 75);
row.(prefix + "_xMedian") = xMedian;
row.(prefix + "_xQ25") = xQ25;
row.(prefix + "_xQ75") = xQ75;
row.(prefix + "_xErrLow") = xMedian-xQ25;
row.(prefix + "_xErrHigh") = xQ75-xMedian;
row.(prefix + "_yMedian") = yMedian;
row.(prefix + "_yQ25") = yQ25;
row.(prefix + "_yQ75") = yQ75;
row.(prefix + "_yErrLow") = yMedian-yQ25;
row.(prefix + "_yErrHigh") = yQ75-yMedian;
end

function validate_wavelength_filter(source)
if ~isfield(source, "wavelengthFilter")
    error("Figure09Depth:MissingFilterProvenance", ...
        "Figure 8 data do not record the real-case wavelength filter.");
end
filter = source.wavelengthFilter;
required = ["minWavelength_m", "maxWavelength_m", ...
    "minAutocorrPeakR", "madScaleFactor", ...
    "scaledMadMultiplier", "minProfilesForMad"];
if ~isstruct(filter) || any(~isfield(filter, cellstr(required))) || ...
        filter.minWavelength_m ~= 0.5 || filter.maxWavelength_m ~= 7 || ...
        filter.minAutocorrPeakR ~= 0.10 || ...
        filter.madScaleFactor ~= 1.4826 || ...
        filter.scaledMadMultiplier ~= 3.5 || ...
        filter.minProfilesForMad ~= 8
    error("Figure09Depth:UnexpectedWavelengthFilter", ...
        "Figure 8 data do not use the agreed real-case wavelength filter.");
end
end
