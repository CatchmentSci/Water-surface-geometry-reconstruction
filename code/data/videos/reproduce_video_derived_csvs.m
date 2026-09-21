function outputs = reproduce_video_derived_csvs(archiveRoot, outputRoot)
%REPRODUCE_VIDEO_DERIVED_CSVS Rebuild canonical real-video CSV products.
%   OUTPUTS = REPRODUCE_VIDEO_DERIVED_CSVS(ARCHIVEROOT, OUTPUTROOT)
%   recalculates the selected-map profiles and projected KLT-IV velocities
%   from deposited checkpoints and inputs. It writes the canonical CSVs
%   below OUTPUTROOT using the same videos/ hierarchy as the archive.
arguments
    archiveRoot (1,1) string
    outputRoot (1,1) string = ""
end

archiveRoot = canonical_path(archiveRoot);
if strlength(outputRoot) == 0
    outputRoot = fullfile(pwd, 'reproduced_video_csvs');
end
outputRoot = canonical_path(outputRoot);

inputFolder = fullfile(archiveRoot, 'videos', 'inputs');
lookupFile = fullfile(inputFolder, 'klt_analysis_case_lookup.tsv');
if ~isfile(lookupFile)
    error('Case lookup was not found: %s', lookupFile);
end

lookup = readtable(lookupFile, 'FileType', 'text', 'Delimiter', '\t', ...
    'CommentStyle', '#', 'TextType', 'string', ...
    'VariableNamingRule', 'preserve');
lookup.filename_in = strtrim(string(lookup.filename_in));
if isnumeric(lookup.sweep_value)
    lookup.sweep_value = double(lookup.sweep_value);
else
    lookup.sweep_value = str2double(string(lookup.sweep_value));
end
lookup = sortrows(lookup, 'sweep_value', 'descend');
lookup = lookup(lookup.sweep_value < 124, :);
if height(lookup) ~= 13
    error('Expected 13 publication cases; found %d.', height(lookup));
end
lookup.case = "R" + string((1:13).');

outputInputs = fullfile(outputRoot, 'videos', 'inputs');
outputProfiles = fullfile(outputRoot, 'videos', 'derived', 'profiles');
outputSummaries = fullfile(outputRoot, 'videos', 'derived', 'summaries');
make_folder(outputInputs);
make_folder(outputProfiles);
make_folder(outputSummaries);

temporaryAuxiliary = string(tempname);
mkdir(temporaryAuxiliary);
cleanup = onCleanup(@() remove_temporary_folder(temporaryAuxiliary)); %#ok<NASGU>

caseOutputs = cell(height(lookup), 1);
velocityRows = cell(height(lookup), 1);
for ii = 1:height(lookup)
    baseName = lookup.filename_in(ii);
    fprintf('\n[%d/%d] Rebuilding %s (%s)\n', ...
        ii, height(lookup), lookup.case(ii), baseName);
    caseOutputs{ii} = generate_case_transect_analysis( ...
        archiveRoot, baseName, temporaryAuxiliary);
    T = caseOutputs{ii}.transectAnalysisTable;
    n = height(T);
    velocityRows{ii} = table( ...
        repmat(lookup.case(ii), n, 1), ...
        double(T.profileIndex), double(T.rowCoord), ...
        double(T.initialVelocityTrackedMedian_mps), ...
        double(T.acceptedVelocityTrackedMedian_mps), ...
        'VariableNames', {'case','profileIndex','rowCoord', ...
        'U_initial','U_accepted'});
end

perTransectTable = vertcat(velocityRows{:});
perTransectFile = fullfile(outputInputs, ...
    'per_transect_initial_accepted.csv');
writetable(perTransectTable, perTransectFile);

profileOutputs = export_selected_profile_csvs( ...
    archiveRoot, temporaryAuxiliary, outputProfiles);
acceptedSummaryFile = export_real_accepted_wsg_summary( ...
    outputProfiles, outputSummaries);

outputs = struct;
outputs.outputRoot = char(outputRoot);
outputs.perTransectFile = char(perTransectFile);
outputs.profileFolder = char(outputProfiles);
outputs.acceptedSummaryFile = char(acceptedSummaryFile);
outputs.nCases = height(lookup);
outputs.nVelocityRows = height(perTransectTable);
outputs.nProfileFiles = profileOutputs.nFilesWritten;

fprintf('\nCanonical CSV rebuild complete:\n');
fprintf('  %s\n', perTransectFile);
fprintf('  %s (%d files)\n', outputProfiles, outputs.nProfileFiles);
fprintf('  %s\n', acceptedSummaryFile);
end

function p = canonical_path(p)
p = string(java.io.File(p).getCanonicalPath());
end

function make_folder(folder)
if ~isfolder(folder)
    mkdir(folder);
end
end

function remove_temporary_folder(folder)
if isfolder(folder)
    rmdir(folder, 's');
end
end
