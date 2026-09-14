function outputs = generate_table_c1(archiveRoot, outputFolder)
%GENERATE_TABLE_C1 Reproduce the real-case wave and hydraulic-results table.
%   OUTPUTS = GENERATE_TABLE_C1(ARCHIVEROOT, OUTPUTFOLDER) combines the
%   deposited accepted-autocorrelation velocity summary, accepted-map
%   amplitudes, and inverse-depth summary for the stable R1-R11 cases.

arguments
    archiveRoot (1, 1) string
    outputFolder (1, 1) string = fullfile(fileparts(mfilename('fullpath')), "output")
end

outputRoot = fullfile(string(archiveRoot), "videos", "outputs");
velocityFile = fullfile(outputRoot, ...
    "wse_autocorrelation_velocity_method_sensitivity_summary.mat");
amplitudeFile = fullfile(outputRoot, "real_wave_wse_pmusic_summary.csv");
depthFile = fullfile(outputRoot, ...
    "wse_autocorrelation_uniform_linear_power_depth_summary.csv");
requiredFiles = [velocityFile, amplitudeFile, depthFile];
for ii = 1:numel(requiredFiles)
    if ~isfile(requiredFiles(ii))
        error("TableC1:MissingInput", ...
            "Required archive file was not found:\n%s", requiredFiles(ii));
    end
end

loaded = load(velocityFile, "autocorrCaseSummaryTable", ...
    "autocorrCaseTransectTables", "autocorrPanelIncluded");
requiredVariables = ["autocorrCaseSummaryTable", ...
    "autocorrCaseTransectTables", "autocorrPanelIncluded"];
missingVariables = setdiff(requiredVariables, string(fieldnames(loaded)));
if ~isempty(missingVariables)
    error("TableC1:InvalidInput", "Velocity MAT-file is missing: %s", ...
        strjoin(missingVariables, ", "));
end

velocity = loaded.autocorrCaseSummaryTable;
transects = loaded.autocorrCaseTransectTables;
included = logical(loaded.autocorrPanelIncluded(:));
amplitude = readtable(amplitudeFile, "VariableNamingRule", "preserve");
depth = readtable(depthFile, "VariableNamingRule", "preserve");

assert_fields(velocity, ["caseNumber", "caseLabel", "deep_yMedian", ...
    "constant_yMedian", "linear_yMedian", "power_yMedian"], ...
    "velocity summary");
assert_fields(amplitude, ["caseNumber", "caseLabel", "amplitudeEst_m"], ...
    "amplitude summary");
assert_fields(depth, ["caseNumber", "caseLabel", "includedInFigure", ...
    "uniform_deep_yMedian", "linear_deep_yMedian", ...
    "power_deep_yMedian"], "depth summary");

if height(velocity) ~= 13 || numel(transects) ~= 13 || ...
        numel(included) ~= 13 || height(amplitude) ~= 13 || height(depth) ~= 13
    error("TableC1:UnexpectedCases", ...
        "Expected 13 aligned real cases in every deposited source.");
end
labels = string(velocity.caseLabel);
if ~isequal(labels, string(amplitude.caseLabel), string(depth.caseLabel)) || ...
        ~isequal(double(velocity.caseNumber), double(amplitude.caseNumber), ...
        double(depth.caseNumber))
    error("TableC1:CaseAlignment", "Deposited source tables are not case-aligned.");
end
if nnz(included) ~= 11 || any(included(12:13)) || ...
        ~isequal(included, logical(depth.includedInFigure))
    error("TableC1:InclusionMismatch", ...
        "Expected the agreed R1-R11 inclusion mask in both analyses.");
end

lambdaEst_m = nan(height(velocity), 1);
for ii = 1:height(velocity)
    T = transects{ii};
    assert_fields(T, ["caseLabel", "acceptedAutocorrWavelength_m", ...
        "commonValid"], sprintf("%s transect table", labels(ii)));
    if any(string(T.caseLabel) ~= labels(ii))
        error("TableC1:CaseAlignment", ...
            "Transect table %d is not aligned with %s.", ii, labels(ii));
    end
    valid = logical(T.commonValid) & ...
        isfinite(double(T.acceptedAutocorrWavelength_m));
    lambdaEst_m(ii) = median(double(T.acceptedAutocorrWavelength_m(valid)));
end

% Independently confirm that each direct wavelength median reproduces the
% deposited deep-water velocity median under U=sqrt(g*lambda/(2*pi)) at the
% table's three-decimal precision. For an even number of transects, taking
% the median before versus after the nonlinear transform can differ slightly.
UfromLambda_mps = sqrt(9.81 .* lambdaEst_m ./ (2*pi));
if any(abs(UfromLambda_mps - double(velocity.deep_yMedian)) >= 5e-4)
    error("TableC1:VelocityMismatch", ...
        "Direct wavelength medians do not reproduce deposited deep-water velocities.");
end

caseNumber = double(velocity.caseNumber(included));
caseLabel = labels(included);
lambdaEst_m = lambdaEst_m(included);
amplitudeEst_m = double(amplitude.amplitudeEst_m(included));
velocityDeep_mps = double(velocity.deep_yMedian(included));
velocityConstant_mps = double(velocity.constant_yMedian(included));
velocityLinear_mps = double(velocity.linear_yMedian(included));
velocityPower_mps = double(velocity.power_yMedian(included));
depthConstant_m = double(depth.uniform_deep_yMedian(included));
depthLinear_m = double(depth.linear_deep_yMedian(included));
depthPower_m = double(depth.power_deep_yMedian(included));

tableData = table(caseNumber, caseLabel, lambdaEst_m, amplitudeEst_m, ...
    velocityDeep_mps, velocityConstant_mps, velocityLinear_mps, ...
    velocityPower_mps, depthConstant_m, depthLinear_m, depthPower_m);

if ~isfolder(outputFolder)
    mkdir(outputFolder);
end
csvFile = fullfile(outputFolder, "table_c1_real_wave_hydraulic_results.csv");
latexFile = fullfile(outputFolder, "table_c1_real_wave_hydraulic_results.tex");
writetable(tableData, csvFile);
write_latex_table(latexFile, tableData);

outputs = struct("csvFile", csvFile, "latexFile", latexFile, ...
    "table", tableData, "velocityFile", velocityFile, ...
    "amplitudeFile", amplitudeFile, "depthFile", depthFile);
end

function write_latex_table(filename, T)
fid = fopen(filename, "w", "n", "UTF-8");
if fid < 0
    error("TableC1:WriteFailed", "Could not open output file: %s", filename);
end
cleanup = onCleanup(@() fclose(fid));

write_line(fid, "\begin{table}[htbp]");
write_line(fid, "\centering");
write_line(fid, "\begin{tabular}{c c c c c c c c c c}");
write_line(fid, "\toprule");
write_line(fid, "& \multicolumn{2}{c}{Estimated}");
write_line(fid, "& \multicolumn{4}{c}{Derived velocity}");
write_line(fid, "& \multicolumn{3}{c}{Derived depth} \\");
write_line(fid, "\cmidrule(lr){2-3}");
write_line(fid, "\cmidrule(lr){4-7}");
write_line(fid, "\cmidrule(lr){8-10}");
write_line(fid, "Case");
write_line(fid, "& $\lambda_{\mathrm{est}}$");
write_line(fid, "& $A_{\mathrm{est}}$");
write_line(fid, "& $U_{\mathrm{deep}}$");
write_line(fid, "& $U_{\mathrm{constant}}$");
write_line(fid, "& $U_{\mathrm{linear}}$");
write_line(fid, "& $U_{\mathrm{power}}$");
write_line(fid, "& $d_{\mathrm{constant}}$");
write_line(fid, "& $d_{\mathrm{linear}}$");
write_line(fid, "& $d_{\mathrm{power}}$ \\");
write_line(fid, "& {[\(\mathrm{m}\)]}");
write_line(fid, "& {[\(\mathrm{m}\)]}");
write_line(fid, "& {[\(\mathrm{m~s^{-1}}\)]}");
write_line(fid, "& {[\(\mathrm{m~s^{-1}}\)]}");
write_line(fid, "& {[\(\mathrm{m~s^{-1}}\)]}");
write_line(fid, "& {[\(\mathrm{m~s^{-1}}\)]}");
write_line(fid, "& {[\(\mathrm{m}\)]}");
write_line(fid, "& {[\(\mathrm{m}\)]}");
write_line(fid, "& {[\(\mathrm{m}\)]} \\");
write_line(fid, "\midrule");

for ii = 1:height(T)
    fields = [sprintf("R$_{%d}$", T.caseNumber(ii)), ...
        format_value(T.lambdaEst_m(ii)), ...
        format_value(T.amplitudeEst_m(ii)), ...
        format_value(T.velocityDeep_mps(ii)), ...
        format_value(T.velocityConstant_mps(ii)), ...
        format_value(T.velocityLinear_mps(ii)), ...
        format_value(T.velocityPower_mps(ii)), ...
        format_value(T.depthConstant_m(ii)), ...
        format_value(T.depthLinear_m(ii)), ...
        format_value(T.depthPower_m(ii))];
    write_line(fid, strjoin(fields, " & ") + " \\");
end

write_line(fid, "\bottomrule");
write_line(fid, "\end{tabular}");
write_line(fid, "\bigskip");
caption = [ ...
    "\caption{Wave-parameter estimates for the real-world cases, ordered " + ...
    "by decreasing discharge; see Table~\ref{Table:real_char} for " + ...
    "hydraulic descriptors. $\lambda_{\mathrm{est}}$ is the median " + ...
    "autocorrelation wavelength estimate returned for the accepted WSG " + ...
    "map. $A_{\mathrm{est}}$ is the median robust amplitude estimate, " + ...
    "calculated for each detrended transect as half the difference between " + ...
    "its 98\textsuperscript{th} and 2\textsuperscript{nd} percentiles. " + ...
    "The derived velocities use the deep-water approximation (Equation~" + ...
    "\ref{eq:3b}), constant profile (Equation~\ref{eq:wavea}), linear " + ...
    "profile (Equation~\ref{eq:waveb}), and power profile (Equation~" + ...
    "\ref{eq:wavec}). Derived depths use the deep solution branch, or the " + ...
    "single solution where only one admissible root exists. R$_{12}$ and " + ...
    "R$_{13}$ are omitted because the selector did not identify stable " + ...
    "wavelength reconstructions.}"];
write_line(fid, caption);
write_line(fid, "\label{Table:real_res}");
write_line(fid, "\end{table}");
end

function text = format_value(value)
if isfinite(value)
    text = string(sprintf("%.3f", value));
else
    text = "--";
end
end

function assert_fields(T, requiredFields, description)
missingFields = setdiff(string(requiredFields), ...
    string(T.Properties.VariableNames));
if ~isempty(missingFields)
    error("TableC1:InvalidInput", "%s is missing field(s): %s", ...
        description, strjoin(missingFields, ", "));
end
end

function write_line(fid, text)
fprintf(fid, "%s\n", text);
end
