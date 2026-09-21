function report = validate_synthetic_summary(referenceRoot, candidateRoot)
%VALIDATE_SYNTHETIC_SUMMARY Check every field used by Figure 5/Table B1.

arguments
    referenceRoot (1, 1) string
    candidateRoot (1, 1) string
end

name = "synthetic_wave_range_autocorr_summary.csv";
referenceFile = fullfile(referenceRoot, "syn", "outputs", name);
candidateFile = fullfile(candidateRoot, name);
A = readtable(referenceFile, "VariableNamingRule", "preserve", "TextType", "string");
B = readtable(candidateFile, "VariableNamingRule", "preserve", "TextType", "string");

fields = ["h_m", "lambda_m", "lambdaEstB_rangeWseAutocorr_m", ...
    "velocityDiffB_rangeWseAutocorr_pct", "A_m", "Aest_m", ...
    "selectedMapIndex", "selectedMapRow", "stopMethod", ...
    "selectionReason", "analysisStatus", "nRangeWseAutocorrProfiles", ...
    "UtrueDeepWater_mps", "UestBAutocorrDeepWater_mps", ...
    "suppressWaveAmpOutputs"];
assert_fields(A, fields);
assert_fields(B, fields);
if height(A) ~= height(B)
    error("SyntheticValidation:RowCount", ...
        "Reference has %d rows; candidate has %d.", height(A), height(B));
end

for field = fields
    a = A.(field);
    b = B.(field);
    if isnumeric(a) || islogical(a)
        if ~isequaln(a, b)
            error("SyntheticValidation:Mismatch", ...
                "Numerical field differs: %s", field);
        end
    elseif ~isequaln(string(a), string(b))
        error("SyntheticValidation:Mismatch", "Text field differs: %s", field);
    end
end

report = table(name, height(A), numel(fields), true, ...
    VariableNames=["file", "rows", "checkedFields", "identical"]);
disp(report);
end

function assert_fields(T, fields)
missing = setdiff(fields, string(T.Properties.VariableNames));
if ~isempty(missing)
    error("SyntheticValidation:MissingField", ...
        "Missing field(s): %s", strjoin(missing, ", "));
end
end
