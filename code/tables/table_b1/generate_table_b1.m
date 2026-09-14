function outputs = generate_table_b1(archiveRoot, outputFolder)
%GENERATE_TABLE_B1 Reproduce the synthetic-case characteristics table.
%   OUTPUTS = GENERATE_TABLE_B1(ARCHIVEROOT, OUTPUTFOLDER) reads the
%   deposited accepted-WSG autocorrelation summary and writes a clean,
%   machine-readable CSV plus the publication-ready Appendix Table B1.

arguments
    archiveRoot (1, 1) string
    outputFolder (1, 1) string = fullfile(fileparts(mfilename('fullpath')), "output")
end

sourceFile = fullfile(string(archiveRoot), "syn", "outputs", ...
    "synthetic_wave_range_autocorr_summary.csv");
if ~isfile(sourceFile)
    error("TableB1:MissingInput", "Required archive file was not found:\n%s", sourceFile);
end

source = readtable(sourceFile, "VariableNamingRule", "preserve");
requiredFields = ["h_m", "lambda_m", ...
    "lambdaEstB_rangeWseAutocorr_m", ...
    "velocityDiffB_rangeWseAutocorr_pct", "A_m", "Aest_m", ...
    "suppressWaveAmpOutputs", "analysisStatus"];
assert_fields(source, requiredFields);
if height(source) ~= 39
    error("TableB1:UnexpectedCases", "Expected 39 cases; found %d.", height(source));
end
if any(string(source.analysisStatus) ~= "ok")
    error("TableB1:InvalidInput", "One or more synthetic cases did not complete successfully.");
end

caseNumber = (1:height(source)).';
caseLabel = "S" + string(caseNumber);
hydraulicDepth_D_m = double(source.h_m);
lambda_m = double(source.lambda_m);
lambdaEst_m = double(source.lambdaEstB_rangeWseAutocorr_m);
deltaU_percent = double(source.velocityDiffB_rangeWseAutocorr_pct);
A_m = double(source.A_m);
Aest_m = double(source.Aest_m);
suppressEstimates = logical(source.suppressWaveAmpOutputs);

% S19 is the single explicitly suppressed case; all other displayed
% estimates must be finite. Preserve genuine high-error estimates.
if nnz(suppressEstimates) ~= 1 || ~suppressEstimates(19)
    error("TableB1:SuppressionMismatch", ...
        "Expected S19 to be the only suppressed synthetic case.");
end
finiteEstimate = isfinite(lambdaEst_m) & isfinite(deltaU_percent) & isfinite(Aest_m);
if any(~finiteEstimate & ~suppressEstimates) || any(finiteEstimate & suppressEstimates)
    error("TableB1:EstimateMismatch", ...
        "Finite estimates do not agree with the deposited suppression flag.");
end

tableData = table(caseNumber, caseLabel, hydraulicDepth_D_m, lambda_m, lambdaEst_m, ...
    deltaU_percent, A_m, Aest_m, suppressEstimates, ...
    VariableNames=["caseNumber", "caseLabel", "hydraulicDepth_D_m", "lambda_m", ...
    "lambdaEst_m", "deltaU_percent", "A_m", "Aest_m", ...
    "suppressEstimates"]);

if ~isfolder(outputFolder)
    mkdir(outputFolder);
end
csvFile = fullfile(outputFolder, "table_b1_synthetic_characteristics.csv");
latexFile = fullfile(outputFolder, "table_b1_synthetic_characteristics.tex");
writetable(tableData, csvFile);
write_latex_table(latexFile, tableData);

outputs = struct("csvFile", csvFile, "latexFile", latexFile, ...
    "table", tableData, "sourceFile", sourceFile);
end

function write_latex_table(filename, T)
fid = fopen(filename, "w");
if fid < 0
    error("TableB1:WriteFailed", "Could not open output file: %s", filename);
end
cleanup = onCleanup(@() fclose(fid));

write_line(fid, "\setcounter{table}{0}");
write_line(fid, "\renewcommand{\thetable}{B\arabic{table}}");
write_line(fid, "\begingroup");
write_line(fid, "\centering");
write_line(fid, "\captionsetup{");
write_line(fid, "font=footnotesize,");
write_line(fid, "skip=4pt");
write_line(fid, "}");
caption = [ ...
    "\captionof{table}{Characteristics associated with each numerical " + ...
    "simulation, in which hydraulic depth $D$, wavelength $\lambda$, and " + ...
    "amplitude $A$ were varied. $\lambda_{\mathrm{est}}$ is the median " + ...
    "autocorrelation wavelength estimated from the accepted WSG map. " + ...
    "$\Delta U$ is the signed percentage difference between the " + ...
    "deep-water gravity-wave velocity calculated from the estimated " + ...
    "wavelength and that expected from the prescribed wavelength, using " + ...
    "$U_s=\sqrt{g\lambda/(2\pi)}$. Estimated values are medians of " + ...
    "streamwise transects spaced at 0.5~m in the cross-stream direction.}"];
write_line(fid, caption);
write_line(fid, "\label{Table:syn_res}");
write_line(fid, "\footnotesize");
write_line(fid, "\setlength{\tabcolsep}{4pt}");
write_line(fid, "\renewcommand{\arraystretch}{0.60}");
write_line(fid, "\begin{adjustbox}{");
write_line(fid, "max width=\textwidth,");
write_line(fid, "max totalheight=0.73\textheight,");
write_line(fid, "keepaspectratio,");
write_line(fid, "center");
write_line(fid, "}");
write_line(fid, "\begin{tabular}{@{}c c c c c c c@{}}");
write_line(fid, "\toprule");
write_line(fid, "Case & $D$ & $\lambda$ & $\lambda_{\mathrm{est}}$ & $\Delta U$ & $A$ & $A_{\mathrm{est}}$ \\");
write_line(fid, " & {[\(\mathrm{m}\)]} & {[\(\mathrm{m}\)]} & {[\(\mathrm{m}\)]} & {[\(\%\)]} & {[\(\mathrm{m}\)]} & {[\(\mathrm{m}\)]} \\");
write_line(fid, "\midrule");

for ii = 1:height(T)
    if T.suppressEstimates(ii)
        lambdaEstText = "---";
        deltaUText = "---";
        AestText = "---";
    else
        lambdaEstText = sprintf("%.3f", T.lambdaEst_m(ii));
        deltaUText = sprintf("%+.2f", T.deltaU_percent(ii));
        AestText = sprintf("%.3f", T.Aest_m(ii));
    end
    row = sprintf("S$_{%d}$ & %.2f & %.3f & %s & %s & %.3f & %s \\\\", ...
        T.caseNumber(ii), T.hydraulicDepth_D_m(ii), T.lambda_m(ii), lambdaEstText, ...
        deltaUText, T.A_m(ii), AestText);
    write_line(fid, row);
end

write_line(fid, "\bottomrule");
write_line(fid, "\end{tabular}");
write_line(fid, "\end{adjustbox}");
write_line(fid, "\endgroup");
end

function assert_fields(T, requiredFields)
missingFields = setdiff(requiredFields, string(T.Properties.VariableNames));
if ~isempty(missingFields)
    error("TableB1:InvalidInput", "Input table is missing field(s): %s", ...
        strjoin(missingFields, ", "));
end
end

function write_line(fid, text)
fprintf(fid, "%s\n", text);
end
