function outputFile = export_real_accepted_wsg_summary(profileFolder, outputFolder)
%EXPORT_REAL_ACCEPTED_WSG_SUMMARY Rebuild the compact 13-case WSG summary.
arguments
    profileFolder (1,1) string
    outputFolder (1,1) string
end

files = dir(fullfile(profileFolder, ...
    '*_selected_map_profiles_for_real_batch_summary.csv'));
if numel(files) ~= 13
    error('Expected 13 selected-profile CSV files in %s; found %d.', ...
        profileFolder, numel(files));
end

rows = repmat(empty_row(), numel(files), 1);
for ii = 1:numel(files)
    sourceFile = fullfile(files(ii).folder, files(ii).name);
    T = readtable(sourceFile, 'TextType', 'string', ...
        'VariableNamingRule', 'preserve');
    required = ["caseNumber","caseLabel","filenameLookup", ...
        "sweepValue_Q","checkpointFileName", ...
        "summaryAmplitudeMedian_m","summaryAmplitudeQ25_m", ...
        "summaryAmplitudeQ75_m","summaryNAmplitudeProfiles", ...
        "selectedMapIndex","selectedMapRow","stopMethod", ...
        "selectionReason","fallbackSelection"];
    missing = setdiff(required, string(T.Properties.VariableNames));
    if ~isempty(missing)
        error('%s is missing column(s): %s', sourceFile, ...
            char(strjoin(missing, ', ')));
    end
    assert_constant_columns(T, required, sourceFile);

    R = empty_row();
    R.caseNumber = double(T.caseNumber(1));
    R.caseLabel = string(T.caseLabel(1));
    R.filenameLookup = string(T.filenameLookup(1));
    R.sweepValue_Q = double(T.sweepValue_Q(1));
    R.checkpointFileName = string(T.checkpointFileName(1));
    R.fileFound = true;
    R.amplitudeEst_m = double(T.summaryAmplitudeMedian_m(1));
    R.amplitudeEst_q25_m = double(T.summaryAmplitudeQ25_m(1));
    R.amplitudeEst_q75_m = double(T.summaryAmplitudeQ75_m(1));
    R.nAmplitudeProfiles = double(T.summaryNAmplitudeProfiles(1));
    R.selectedMapIndex = double(T.selectedMapIndex(1));
    R.selectedMapRow = double(T.selectedMapRow(1));
    R.stopMethod = string(T.stopMethod(1));
    R.selectionReason = string(T.selectionReason(1));
    R.fallbackSelection = logical(T.fallbackSelection(1));
    R.analysisStatus = "ok";
    R.failureMessage = "";
    rows(ii) = R;
end

summaryTable = sortrows(struct2table(rows), 'caseNumber');
if ~isequal(double(summaryTable.caseNumber), (1:13).') || ...
        ~isequal(string(summaryTable.caseLabel), "R" + string((1:13).'))
    error('Selected-profile files are not the complete ordered R1-R13 set.');
end
if ~isfolder(outputFolder)
    mkdir(outputFolder);
end
outputFile = fullfile(outputFolder, 'real_accepted_wsg_summary.csv');
writetable(summaryTable, outputFile);
outputFile = char(outputFile);
end

function R = empty_row()
R = struct( ...
    'caseNumber', NaN, ...
    'caseLabel', "", ...
    'filenameLookup', "", ...
    'sweepValue_Q', NaN, ...
    'checkpointFileName', "", ...
    'fileFound', false, ...
    'amplitudeEst_m', NaN, ...
    'amplitudeEst_q25_m', NaN, ...
    'amplitudeEst_q75_m', NaN, ...
    'nAmplitudeProfiles', NaN, ...
    'selectedMapIndex', NaN, ...
    'selectedMapRow', NaN, ...
    'stopMethod', "", ...
    'selectionReason', "", ...
    'fallbackSelection', false, ...
    'analysisStatus', "", ...
    'failureMessage', "");
end

function assert_constant_columns(T, names, sourceFile)
for jj = 1:numel(names)
    value = T.(names(jj));
    if isnumeric(value) || islogical(value)
        same = all(arrayfun(@(kk) isequaln(value(1), value(kk)), ...
            2:numel(value)));
    else
        same = numel(unique(string(value))) == 1;
    end
    if ~same
        error('%s has non-constant case field %s.', sourceFile, names(jj));
    end
end
end
