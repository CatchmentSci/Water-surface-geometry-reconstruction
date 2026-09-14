function outputs = generate_table_03(archiveRoot, outputFolder)
%GENERATE_TABLE_03 Reproduce the real-case flow-characteristics table.
%   OUTPUTS = GENERATE_TABLE_03(ARCHIVEROOT, OUTPUTFOLDER) combines the
%   deposited hydraulic statistics with seeding densities recalculated
%   from the 13 real-case solver checkpoints. It writes a machine-readable
%   CSV and publication-ready LaTeX table.

arguments
    archiveRoot (1, 1) string
    outputFolder (1, 1) string = fullfile(fileparts(mfilename('fullpath')), "output")
end

archiveRoot = string(archiveRoot);
outputFolder = string(outputFolder);
inputFolder = fullfile(archiveRoot, "videos", "inputs");
checkpointFolder = fullfile(archiveRoot, "videos", "outputs");

hydraulicFile = fullfile(inputFolder, "Dart_video_hydraulic_statistics.csv");
lookupFile = fullfile(inputFolder, "klt_analysis_case_lookup.tsv");
sweepFile = fullfile(inputFolder, "sweep_limits.csv");
selectionFile = fullfile(checkpointFolder, "real_wave_wse_pmusic_summary.csv");
requiredFiles = [hydraulicFile, lookupFile, sweepFile, selectionFile];
for ii = 1:numel(requiredFiles)
    if ~isfile(requiredFiles(ii))
        error("Table03:MissingInput", ...
            "Required archive file was not found:\n%s", requiredFiles(ii));
    end
end

hydraulic = readtable(hydraulicFile, "VariableNamingRule", "preserve");
requiredHydraulicFields = ["FolderName", "Discharge_m3s", ...
    "FlowExceedance_percent", "MeanVelocity_mps", "AverageDepth_m", ...
    "FroudeNumber", "ReynoldsNumber_Rh", "NearestRecordOffset_minutes"];
assert_fields(hydraulic, requiredHydraulicFields, "hydraulic table");

% The paper's R1-R13 set comprises the 13 analysed cases below 124 m3 s-1,
% ordered by decreasing discharge.
hydraulic = hydraulic(double(hydraulic.Discharge_m3s) < 124, :);
hydraulic = sortrows(hydraulic, "Discharge_m3s", "descend");
if height(hydraulic) ~= 13
    error("Table03:UnexpectedCases", ...
        "Expected 13 hydraulic cases below 124 m3 s-1; found %d.", ...
        height(hydraulic));
end

lookup = read_lookup_table(lookupFile);
assert_fields(lookup, ["filename_in", "sweep_value"], "case lookup");
sweepLimits = readtable(sweepFile, "VariableNamingRule", "preserve");
assert_fields(sweepLimits, ["videoNumber", "row_near", "row_far"], ...
    "sweep-limits table");
selection = readtable(selectionFile, "VariableNamingRule", "preserve");
assert_fields(selection, ["caseLabel", "filenameLookup", ...
    "selectedMapIndex"], "accepted-map summary");

nCases = height(hydraulic);
caseNumber = (1:nCases).';
caseLabel = "R" + string(caseNumber);
sweepValue = round(double(hydraulic.Discharge_m3s));
videoIdentifier = strings(nCases, 1);
meanDensity_pts_m2 = nan(nCases, 1);
densityDetails = repmat(empty_density_details(), nCases, 1);

for ii = 1:nCases
    lookupRows = find(double(lookup.sweep_value) == sweepValue(ii));
    if numel(lookupRows) ~= 1
        error("Table03:CaseLookup", ...
            "Expected one lookup row for sweep value %g; found %d.", ...
            sweepValue(ii), numel(lookupRows));
    end
    videoIdentifier(ii) = string(lookup.filename_in(lookupRows));
    checkpointFile = fullfile(checkpointFolder, ...
        videoIdentifier(ii) + "_checkpoint.mat");
    if ~isfile(checkpointFile)
        error("Table03:MissingCheckpoint", ...
            "Required checkpoint was not found:\n%s", checkpointFile);
    end

    limitRows = find(double(sweepLimits.videoNumber) == sweepValue(ii));
    if numel(limitRows) ~= 1
        error("Table03:SweepLimits", ...
            "Expected one sweep-limits row for value %g; found %d.", ...
            sweepValue(ii), numel(limitRows));
    end
    rowLimits = sort([double(sweepLimits.row_near(limitRows)), ...
        double(sweepLimits.row_far(limitRows))]);

    selectionRows = find(string(selection.caseLabel) == caseLabel(ii));
    if numel(selectionRows) ~= 1 || ...
            string(selection.filenameLookup(selectionRows)) ~= videoIdentifier(ii)
        error("Table03:MapSelection", ...
            "Accepted-map summary does not uniquely match %s.", caseLabel(ii));
    end
    selectedMapIndex = double(selection.selectedMapIndex(selectionRows));

    loaded = load(checkpointFile, "checkpoint");
    if ~isfield(loaded, "checkpoint") || ~isstruct(loaded.checkpoint)
        error("Table03:InvalidCheckpoint", ...
            "Checkpoint variable was not found in:\n%s", checkpointFile);
    end
    details = calculate_seeding_density( ...
        loaded.checkpoint, rowLimits, selectedMapIndex);
    details.checkpointFile = checkpointFile;
    details.caseLabel = caseLabel(ii);
    details.sweepValue = sweepValue(ii);
    densityDetails(ii) = details;
    meanDensity_pts_m2(ii) = details.meanDensity_pts_m2;
    fprintf("%s: %.6f pts m^-2\n", caseLabel(ii), meanDensity_pts_m2(ii));
end

% Use conventional nearest-integer rounding for the displayed table while
% preserving full-precision density values in the CSV output.
densityDisplay_pts_m2 = round(meanDensity_pts_m2);
maximumNearestOffset_minutes = max( ...
    double(hydraulic.NearestRecordOffset_minutes), [], "omitnan");

tableData = table(caseNumber, caseLabel, videoIdentifier, sweepValue, ...
    double(hydraulic.Discharge_m3s), ...
    double(hydraulic.FlowExceedance_percent), ...
    double(hydraulic.AverageDepth_m), ...
    double(hydraulic.MeanVelocity_mps), ...
    double(hydraulic.FroudeNumber), ...
    double(hydraulic.ReynoldsNumber_Rh), ...
    meanDensity_pts_m2, densityDisplay_pts_m2, ...
    VariableNames=["caseNumber", "caseLabel", "videoIdentifier", ...
    "sweepValue", "discharge_Q_m3ps", "flowExceedance_percent", ...
    "hydraulicDepth_D_m", "sectionMeanVelocity_u_mps", ...
    "froudeNumber", "reynoldsNumber_Rh", ...
    "meanSeedingDensity_pts_m2", "displaySeedingDensity_pts_m2"]);

if ~isfolder(outputFolder)
    mkdir(outputFolder);
end
csvFile = fullfile(outputFolder, "table_03_flow_characteristics.csv");
latexFile = fullfile(outputFolder, "table_03_flow_characteristics.tex");
writetable(tableData, csvFile);
write_latex_table(latexFile, tableData, maximumNearestOffset_minutes);

outputs = struct;
outputs.csvFile = csvFile;
outputs.latexFile = latexFile;
outputs.table = tableData;
outputs.densityDetails = densityDetails;
outputs.maximumNearestRecordOffset_minutes = maximumNearestOffset_minutes;
end

function lookup = read_lookup_table(filename)
lines = readlines(filename);
lines = lines(strlength(strtrim(lines)) > 0);
lines = lines(~startsWith(strtrim(lines), "#"));
temporaryFile = string(tempname) + ".tsv";
cleanup = onCleanup(@() delete_if_present(temporaryFile));
writelines(lines, temporaryFile);
lookup = readtable(temporaryFile, "FileType", "text", ...
    "Delimiter", "\t", "VariableNamingRule", "preserve");
end

function details = calculate_seeding_density(checkpoint, rowLimits, mapIndex)
aa = 1;
if ~isfield(checkpoint, "wse_map") || ~iscell(checkpoint.wse_map)
    error("Table03:InvalidCheckpoint", "checkpoint.wse_map is missing.");
end
mapIndex = round(mapIndex);
if ~isfinite(mapIndex) || mapIndex < 1 || ...
        mapIndex > size(checkpoint.wse_map, 2) || ...
        isempty(checkpoint.wse_map{aa, mapIndex})
    error("Table03:InvalidCheckpoint", ...
        "Accepted WSG map index %g is unavailable.", mapIndex);
end
W = checkpoint.wse_map{aa, mapIndex};
[nRows, nColumns] = size(W);
if rowLimits(1) < 1 || rowLimits(2) > nRows || ...
        any(rowLimits ~= round(rowLimits))
    error("Table03:InvalidRowLimits", ...
        "Selected rows %g:%g are outside the WSG map.", ...
        rowLimits(1), rowLimits(2));
end

if isfield(checkpoint, "dx") && isfinite(checkpoint.dx)
    dx = double(checkpoint.dx);
else
    dx = median(abs(diff(double(checkpoint.xi))), "omitnan");
end
if isfield(checkpoint, "dy") && isfinite(checkpoint.dy)
    dy = double(checkpoint.dy);
else
    dy = median(abs(diff(double(checkpoint.yi))), "omitnan");
end
cellArea_m2 = abs(dx * dy);
if ~isfinite(cellArea_m2) || cellArea_m2 <= 0
    error("Table03:InvalidCheckpoint", ...
        "A positive WSG cell area could not be determined.");
end

if isfield(checkpoint, "dzMapHist") && ~isempty(checkpoint.dzMapHist)
    dzHistory = checkpoint.dzMapHist;
elseif isfield(checkpoint, "metrics_static") && ...
        isfield(checkpoint.metrics_static, "dzMapHist") && ...
        ~isempty(checkpoint.metrics_static.dzMapHist)
    dzHistory = checkpoint.metrics_static.dzMapHist;
else
    error("Table03:InvalidCheckpoint", ...
        "No dzMapHist was found to define the adjusted-cell mask.");
end

maxIteration = min(numel(dzHistory), max(mapIndex - 1, 0));
adjustedMask = false(size(W));
for kk = 1:maxIteration
    dz = dzHistory{kk};
    if ~isempty(dz) && isequal(size(dz), size(W))
        adjustedMask = adjustedMask | (isfinite(dz) & abs(dz) > 0);
    end
end
adjustedMask = bwareaopen(adjustedMask, 20);
if any(adjustedMask(:))
    adjustedMask = bwareafilt(adjustedMask, 1);
end
rowMask = false(size(W));
rowMask(rowLimits(1):rowLimits(2), :) = true;
adjustedMask = adjustedMask & rowMask;
if ~any(adjustedMask(:))
    error("Table03:InvalidCheckpoint", ...
        "Adjusted-cell mask is empty after applying the case row limits.");
end

countHistory = {};
if isfield(checkpoint, "metrics_static") && ...
        isfield(checkpoint.metrics_static, "debug") && ...
        isfield(checkpoint.metrics_static.debug, "candidateCountAStartHist")
    countHistory = checkpoint.metrics_static.debug.candidateCountAStartHist;
elseif isfield(checkpoint, "metrics") && ...
        isfield(checkpoint.metrics, "debug") && ...
        isfield(checkpoint.metrics.debug, "candidateCountAStartHist")
    countHistory = checkpoint.metrics.debug.candidateCountAStartHist;
end
if isempty(countHistory)
    error("Table03:InvalidCheckpoint", ...
        "candidateCountAStartHist is missing; raw point density cannot be calculated.");
end

pointCountMap = [];
for kk = 1:min(numel(countHistory), maxIteration)
    candidateMap = countHistory{kk};
    if ~isempty(candidateMap) && isequal(size(candidateMap), size(W))
        pointCountMap = double(candidateMap);
        break
    end
end
if isempty(pointCountMap)
    error("Table03:InvalidCheckpoint", ...
        "No valid candidate-count map was found.");
end

nAdjustedCells = nnz(adjustedMask);
adjustedArea_m2 = nAdjustedCells * cellArea_m2;
totalEquivalentPoints = sum(pointCountMap(adjustedMask), "omitnan");
details = empty_density_details();
details.mapIndex = mapIndex;
details.rowNear = rowLimits(1);
details.rowFar = rowLimits(2);
details.nRows = nRows;
details.nColumns = nColumns;
details.dx_m = dx;
details.dy_m = dy;
details.cellArea_m2 = cellArea_m2;
details.nAdjustedCells = nAdjustedCells;
details.adjustedArea_m2 = adjustedArea_m2;
details.totalEquivalentPoints = totalEquivalentPoints;
details.meanDensity_pts_m2 = totalEquivalentPoints / adjustedArea_m2;
end

function details = empty_density_details()
details = struct("checkpointFile", "", "caseLabel", "", ...
    "sweepValue", NaN, "mapIndex", NaN, "rowNear", NaN, ...
    "rowFar", NaN, "nRows", NaN, "nColumns", NaN, "dx_m", NaN, ...
    "dy_m", NaN, "cellArea_m2", NaN, "nAdjustedCells", NaN, ...
    "adjustedArea_m2", NaN, "totalEquivalentPoints", NaN, ...
    "meanDensity_pts_m2", NaN);
end

function write_latex_table(filename, T, maximumNearestOffset_minutes)
fid = fopen(filename, "w");
if fid < 0
    error("Table03:WriteFailed", "Could not open output file: %s", filename);
end
cleanup = onCleanup(@() fclose(fid));
write_line(fid, "\begin{table}[htbp]");
write_line(fid, "\centering");
write_line(fid, "\begin{tabular}{c c c c c c @{\hspace{0.35cm}} c @{\hspace{0.35cm}} c}");
write_line(fid, "\toprule");
write_line(fid, "Case & $Q$ & Exceedance & $D$ & $u$ & $Fr$ & $Re$ & Density \\");
write_line(fid, " & {[\(\mathrm{m^3~s^{-1}}\)]} & {[\(\%\)]} & {[\(\mathrm{m}\)]} & {[\(\mathrm{m~s^{-1}}\)]} & {[\(-\)]} & {[\(-\)]} & {[\(\mathrm{pts~m^{-2}}\)]} \\");
write_line(fid, "\midrule");
for ii = 1:height(T)
    reText = scientific_latex(T.reynoldsNumber_Rh(ii));
    rowFormat = "R$_{%d}$ & %.2f & %.2f & %.2f & %.2f & %.2f & " + ...
        "\\(%s\\) & %d \\\\";
    row = sprintf(rowFormat, T.caseNumber(ii), T.discharge_Q_m3ps(ii), ...
        T.flowExceedance_percent(ii), T.hydraulicDepth_D_m(ii), ...
        T.sectionMeanVelocity_u_mps(ii), T.froudeNumber(ii), reText, ...
        T.displaySeedingDensity_pts_m2(ii));
    write_line(fid, row);
end
write_line(fid, "\bottomrule");
write_line(fid, "\end{tabular}");
write_line(fid, "\bigskip");
captionFormat = "\\caption{Flow characteristics associated with each " + ...
    "real-world case analysed. Where video acquisition and river-stage " + ...
    "measurements were not synchronous, stage and discharge were linearly " + ...
    "interpolated between adjacent records; the maximum offset to the " + ...
    "nearest supporting record was %.0f~min. Flow-exceedance statistics " + ...
    "were calculated from the long-term record spanning 1958--2026. " + ...
    "$D$ is hydraulic depth, $u$ is section-averaged velocity, and density " + ...
    "is the mean seeding density across the analysed domain.}";
caption = sprintf(captionFormat, ...
    maximumNearestOffset_minutes);
write_line(fid, caption);
write_line(fid, "\label{Table:real_char}");
write_line(fid, "\end{table}");
end

function text = scientific_latex(value)
exponent = floor(log10(abs(value)));
mantissa = value / 10^exponent;
text = sprintf("%.1f \\times 10^{%d}", mantissa, exponent);
end

function assert_fields(T, requiredFields, description)
missingFields = setdiff(requiredFields, string(T.Properties.VariableNames));
if ~isempty(missingFields)
    error("Table03:InvalidInput", "%s is missing field(s): %s", ...
        description, strjoin(missingFields, ", "));
end
end

function write_line(fid, text)
fprintf(fid, "%s\n", text);
end

function delete_if_present(filename)
if isfile(filename)
    delete(filename);
end
end
