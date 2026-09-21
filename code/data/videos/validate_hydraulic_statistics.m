function report = validate_hydraulic_statistics(referenceRoot, candidateRoot)
%VALIDATE_HYDRAULIC_STATISTICS Compare a rebuilt hydraulic table exactly.

arguments
    referenceRoot (1, 1) string
    candidateRoot (1, 1) string
end

name = "Dart_video_hydraulic_statistics.csv";
referenceFile = fullfile(referenceRoot, "videos", "inputs", name);
candidateFile = fullfile(candidateRoot, name);
A = readtable(referenceFile, "VariableNamingRule", "preserve", "TextType", "string");
B = readtable(candidateFile, "VariableNamingRule", "preserve", "TextType", "string");

if ~isequal(A.Properties.VariableNames, B.Properties.VariableNames)
    error("HydraulicValidation:Columns", "Column names or order differ.");
end
if height(A) ~= height(B)
    error("HydraulicValidation:Rows", ...
        "Reference has %d rows; candidate has %d.", height(A), height(B));
end

for ii = 1:width(A)
    name_i = A.Properties.VariableNames{ii};
    a = A.(name_i);
    b = B.(name_i);
    if isnumeric(a) || islogical(a) || isdatetime(a)
        same = isequaln(a, b);
    else
        same = isequaln(string(a), string(b));
    end
    if ~same
        error("HydraulicValidation:Mismatch", "Field differs: %s", name_i);
    end
end

report = table(name, height(A), width(A), true, ...
    VariableNames=["file", "rows", "columns", "identical"]);
disp(report);
end
