function report = validate_video_derived_csvs(referenceRoot, candidateRoot)
%VALIDATE_VIDEO_DERIVED_CSVS Compare rebuilt CSVs with deposited products.
arguments
    referenceRoot (1,1) string
    candidateRoot (1,1) string
end

relativeFiles = [ ...
    "videos/inputs/per_transect_initial_accepted.csv"; ...
    "videos/derived/summaries/real_accepted_wsg_summary.csv"];

referenceProfiles = dir(fullfile(referenceRoot, 'videos', 'derived', ...
    'profiles', '*_selected_map_profiles_for_real_batch_summary.csv'));
candidateProfiles = dir(fullfile(candidateRoot, 'videos', 'derived', ...
    'profiles', '*_selected_map_profiles_for_real_batch_summary.csv'));
referenceNames = sort(string({referenceProfiles.name}).');
candidateNames = sort(string({candidateProfiles.name}).');
if ~isequal(referenceNames, candidateNames) || numel(referenceNames) ~= 13
    error('Reference and candidate profile-file inventories differ.');
end
relativeFiles = [relativeFiles; ...
    "videos/derived/profiles/" + referenceNames];

rows = repmat(struct('file',"",'rows',0,'columns',0,'status',""), ...
    numel(relativeFiles), 1);
for ii = 1:numel(relativeFiles)
    relativeFile = relativeFiles(ii);
    referenceFile = fullfile(referenceRoot, relativeFile);
    candidateFile = fullfile(candidateRoot, relativeFile);
    if ~isfile(referenceFile) || ~isfile(candidateFile)
        error('Missing comparison file: %s', relativeFile);
    end
    A = readtable(referenceFile, 'TextType', 'string', ...
        'VariableNamingRule', 'preserve');
    B = readtable(candidateFile, 'TextType', 'string', ...
        'VariableNamingRule', 'preserve');
    compare_tables(A, B, relativeFile);
    rows(ii).file = relativeFile;
    rows(ii).rows = height(A);
    rows(ii).columns = width(A);
    rows(ii).status = "identical";
end
report = struct2table(rows);
fprintf('Validated %d canonical CSV files with no differences.\n', ...
    height(report));
end

function compare_tables(A, B, label)
if height(A) ~= height(B) || width(A) ~= width(B) || ...
        ~isequal(string(A.Properties.VariableNames), ...
        string(B.Properties.VariableNames))
    error('%s has a different table shape or column inventory.', label);
end

for jj = 1:width(A)
    name = A.Properties.VariableNames{jj};
    a = A.(name);
    b = B.(name);
    if isnumeric(a) && isnumeric(b)
        bothNan = isnan(a) & isnan(b);
        difference = abs(double(a) - double(b));
        scale = max(1, max(abs([double(a(:)); double(b(:))]), [], ...
            'omitnan'));
        equal = bothNan | difference <= 1e-12 .* scale;
    elseif islogical(a) && islogical(b)
        equal = a == b;
    else
        equal = string(a) == string(b);
    end
    if ~all(equal, 'all')
        first = find(~equal, 1, 'first');
        error('%s differs in column %s at row %d.', label, name, first);
    end
end
end
