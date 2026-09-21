function outputs = build_figure_08_summary(archiveRoot, outputFolder)
%BUILD_FIGURE_08_SUMMARY Build Figure 8 data from canonical deposited CSVs.
arguments
    archiveRoot (1,1) string
    outputFolder (1,1) string
end

velocityFile = fullfile(archiveRoot, 'videos', 'inputs', ...
    'per_transect_initial_accepted.csv');
profileFolder = fullfile(archiveRoot, 'videos', 'derived', 'profiles');
if ~isfile(velocityFile) || ~isfolder(profileFolder)
    error('Canonical velocity or profile inputs were not found below %s.', ...
        archiveRoot);
end
if ~isfolder(outputFolder)
    mkdir(outputFolder);
end

velocity = readtable(velocityFile, 'TextType', 'string', ...
    'VariableNamingRule', 'preserve');
requiredVelocity = ["case","profileIndex","U_initial","U_accepted"];
missing = setdiff(requiredVelocity, ...
    string(velocity.Properties.VariableNames));
if ~isempty(missing)
    error('Velocity CSV is missing column(s): %s', ...
        char(strjoin(missing, ', ')));
end

methods = [ ...
    struct('key','constant','label','Constant profile','marker','o'), ...
    struct('key','linear','label','Linear profile','marker','s'), ...
    struct('key','power','label','Power profile','marker','d')];
keys = ["deep","constant","linear","power"];
caseTransects = cell(13, 1);
summaryRows = repmat(empty_case_summary(keys), 13, 1);

profileFiles = dir(fullfile(profileFolder, ...
    '*_selected_map_profiles_for_real_batch_summary.csv'));
if numel(profileFiles) ~= 13
    error('Expected 13 selected-profile CSV files; found %d.', ...
        numel(profileFiles));
end

for ii = 1:numel(profileFiles)
    P = readtable(fullfile(profileFiles(ii).folder, profileFiles(ii).name), ...
        'TextType', 'string', 'VariableNamingRule', 'preserve');
    requiredProfile = ["caseNumber","caseLabel","sweepValue_Q", ...
        "profileIndex","crossSectionDepth_m", ...
        "autocorrWavelength_m","autocorrPeakR"];
    missing = setdiff(requiredProfile, string(P.Properties.VariableNames));
    if ~isempty(missing)
        error('%s is missing column(s): %s', profileFiles(ii).name, ...
            char(strjoin(missing, ', ')));
    end
    caseNumber = double(P.caseNumber(1));
    if ~isscalar(caseNumber) || caseNumber < 1 || caseNumber > 13
        error('%s has an invalid case number.', profileFiles(ii).name);
    end
    P = sortrows(P, 'profileIndex');
    caseLabel = string(P.caseLabel(1));
    V = velocity(string(velocity.case) == caseLabel, :);
    V = sortrows(V, 'profileIndex');
    if ~isequal(double(P.profileIndex), double(V.profileIndex))
        error('Velocity/profile keys do not align for %s.', caseLabel);
    end

    lambda = accepted_autocorr_wavelengths( ...
        double(P.autocorrWavelength_m), double(P.autocorrPeakR));
    depth = double(P.crossSectionDepth_m);
    initialVelocity = double(V.U_initial);
    observed = double(V.U_accepted);
    [Udeep,Uconstant,Ulinear,Upower] = calculate_velocities( ...
        lambda, depth, 9.81, 0.85);
    commonValid = isfinite(observed) & observed > 0 & ...
        isfinite(lambda) & lambda > 0 & isfinite(depth) & depth > 0 & ...
        isfinite(Udeep) & Udeep > 0 & ...
        isfinite(Uconstant) & Uconstant > 0 & ...
        isfinite(Ulinear) & Ulinear > 0 & ...
        isfinite(Upower) & Upower > 0;

    n = height(P);
    T = table(repmat(caseLabel,n,1), ...
        repmat(double(P.sweepValue_Q(1)),n,1), double(P.profileIndex), ...
        observed, depth, lambda, Udeep, Uconstant, Ulinear, Upower, ...
        commonValid, initialVelocity, ...
        'VariableNames', {'caseLabel','discharge_Q_m3ps','profileIndex', ...
        'observedVelocity_mps','crossSectionDepth_m', ...
        'acceptedAutocorrWavelength_m','velocityDeepWater_mps', ...
        'velocityConstantProfile_mps','velocityLinearProfile_mps', ...
        'velocityPowerProfile_mps','commonValid', ...
        'initialVelocityTrackedMedian_mps'});
    caseTransects{caseNumber} = T;

    S = empty_case_summary(keys);
    S.caseNumber = caseNumber;
    S.caseLabel = caseLabel;
    S.discharge_Q_m3ps = double(P.sweepValue_Q(1));
    S.nTransects = n;
    S.nCommonValidTransects = nnz(commonValid);
    estimates = struct('deep',Udeep,'constant',Uconstant, ...
        'linear',Ulinear,'power',Upower);
    for kk = 1:numel(keys)
        key = char(keys(kk));
        estimate = estimates.(key);
        S = add_pair_summary(S, initialVelocity(commonValid), ...
            estimate(commonValid), key);
    end
    summaryRows(caseNumber) = S;
end

if any(cellfun(@isempty, caseTransects))
    error('The profile files did not resolve to a complete R1-R13 set.');
end
autocorrCaseSummaryTable = struct2table(summaryRows);
autocorrCaseTransectTables = caseTransects;
autocorrSensitivitySummaryTable = summarise_sensitivity( ...
    autocorrCaseSummaryTable, autocorrCaseTransectTables, methods);
autocorrPanelIncluded = ~ismember( ...
    string(autocorrCaseSummaryTable.caseLabel), ["R12","R13"]);
autocorrSensitivitySummaryTable.includedInSensitivityPanel = ...
    autocorrPanelIncluded;
autocorrSensitivitySummaryTable.sensitivityExclusionReason = ...
    strings(13,1);
autocorrSensitivitySummaryTable.sensitivityExclusionReason( ...
    ~autocorrPanelIncluded) = "Rejected erroneous wavelength reconstruction";
sensitivityMethods = methods;

summaryMatFile = fullfile(outputFolder, ...
    'wse_autocorrelation_velocity_method_sensitivity_summary.mat');
summaryCsvFile = fullfile(outputFolder, ...
    'wse_autocorrelation_velocity_method_sensitivity_summary.csv');
sourceCsvFiles = struct('velocity',char(velocityFile), ...
    'profiles',char(profileFolder));
save(summaryMatFile, 'autocorrCaseSummaryTable', ...
    'autocorrCaseTransectTables', 'autocorrSensitivitySummaryTable', ...
    'autocorrPanelIncluded', 'sensitivityMethods', 'sourceCsvFiles');
writetable(autocorrSensitivitySummaryTable, summaryCsvFile);

outputs = struct('matFile',char(summaryMatFile), ...
    'csvFile',char(summaryCsvFile), 'nCases',13, ...
    'nTransects',sum(cellfun(@height,caseTransects)));
end

function lambda = accepted_autocorr_wavelengths(lambda, peakR)
lambda = double(lambda(:));
peakR = double(peakR(:));
lambda(~isfinite(lambda) | lambda <= 0 | ~isfinite(peakR) | peakR < 0.10) = NaN;
idx = find(isfinite(lambda));
if numel(idx) >= 8
    values = lambda(idx);
    centre = median(values, 'omitnan');
    scaledMad = 1.4826 .* median(abs(values-centre), 'omitnan');
    if isfinite(scaledMad) && scaledMad > 0
        lambda(idx(abs(values-centre) > 3.5.*scaledMad)) = NaN;
    end
end
end

function [Udeep,Uconstant,Ulinear,Upower] = ...
        calculate_velocities(lambda, depth, g, alpha)
lambda = double(lambda(:)); depth = double(depth(:));
valid = isfinite(lambda) & lambda > 0 & isfinite(depth) & depth > 0;
nValues = numel(lambda);
Udeep=nan(nValues,1); Uconstant=Udeep; Ulinear=Udeep; Upower=Udeep;
k=nan(nValues,1); kh=k;
k(valid)=2*pi./lambda(valid); kh(valid)=k(valid).*depth(valid);
Udeep(valid)=sqrt(g./k(valid));
idx=find(valid);
arg=g.*depth(valid).*tanh(kh(valid))./kh(valid);
ok=isfinite(arg)&arg>0; Uconstant(idx(ok))=sqrt(arg(ok));
m=2.*(1-alpha); den=kh(valid)-m.*tanh(kh(valid));
arg=g.*depth(valid).*tanh(kh(valid))./den;
ok=isfinite(arg)&arg>0&isfinite(den)&den>0;
Ulinear(idx(ok))=sqrt(arg(ok));
n=1./alpha-1; s=sign(0.5-n);
ratio=besseli(s.*(0.5-n),kh(valid),1)./ ...
    besseli(-s.*(0.5+n),kh(valid),1);
arg=g.*depth(valid).*ratio./kh(valid);
ok=isfinite(arg)&isreal(arg)&real(arg)>0;
Upower(idx(ok))=sqrt(real(arg(ok)));
end

function S = empty_case_summary(keys)
S=struct('caseNumber',NaN,'caseLabel',"",'discharge_Q_m3ps',NaN, ...
    'nTransects',0,'nCommonValidTransects',0);
for ii=1:numel(keys)
    key=char(keys(ii));
    names={'nPairs','xMedian','xQ25','xQ75','xErrLow','xErrHigh', ...
        'yMedian','yQ25','yQ75','yErrLow','yErrHigh'};
    for jj=1:numel(names)
        value=NaN; if strcmp(names{jj},'nPairs'), value=0; end
        S.([key '_' names{jj}])=value;
    end
end
end

function S = add_pair_summary(S,x,y,key)
valid=isfinite(x)&isfinite(y); x=x(valid); y=y(valid);
S.([key '_nPairs'])=numel(x); if isempty(x), return; end
xm=median(x); x25=prctile(x,25); x75=prctile(x,75);
ym=median(y); y25=prctile(y,25); y75=prctile(y,75);
S.([key '_xMedian'])=xm; S.([key '_xQ25'])=x25;
S.([key '_xQ75'])=x75; S.([key '_xErrLow'])=xm-x25;
S.([key '_xErrHigh'])=x75-xm; S.([key '_yMedian'])=ym;
S.([key '_yQ25'])=y25; S.([key '_yQ75'])=y75;
S.([key '_yErrLow'])=ym-y25; S.([key '_yErrHigh'])=y75-ym;
end

function tableOut = summarise_sensitivity(caseSummary,caseTables,methods)
rows=repmat(empty_sensitivity(methods),height(caseSummary),1);
fields=struct('constant','velocityConstantProfile_mps', ...
    'linear','velocityLinearProfile_mps','power','velocityPowerProfile_mps');
for ii=1:height(caseSummary)
    T=caseTables{ii}; lambda=T.acceptedAutocorrWavelength_m;
    kh=2*pi.*T.crossSectionDepth_m./lambda; Udeep=T.velocityDeepWater_mps;
    valid=T.commonValid&isfinite(kh)&kh>0&isfinite(Udeep)&Udeep>0;
    S=empty_sensitivity(methods); S.caseNumber=caseSummary.caseNumber(ii);
    S.caseLabel=caseSummary.caseLabel(ii);
    S.discharge_Q_m3ps=caseSummary.discharge_Q_m3ps(ii);
    S.nCommonValidTransects=nnz(valid);
    if any(valid)
        v=kh(valid); med=median(v); q25=prctile(v,25); q75=prctile(v,75);
        S.khMedian=med; S.khQ25=q25; S.khQ75=q75;
        S.khErrLow=med-q25; S.khErrHigh=q75-med;
    end
    for jj=1:numel(methods)
        key=methods(jj).key; U=T.(fields.(key)); ok=valid&isfinite(U)&U>0;
        d=100.*(U(ok)-Udeep(ok))./Udeep(ok);
        S.([key '_nPairs'])=nnz(ok);
        if ~isempty(d)
            med=median(d); q25=prctile(d,25); q75=prctile(d,75);
            S.([key '_differenceMedian_pct'])=med;
            S.([key '_differenceQ25_pct'])=q25;
            S.([key '_differenceQ75_pct'])=q75;
            S.([key '_differenceErrLow_pct'])=med-q25;
            S.([key '_differenceErrHigh_pct'])=q75-med;
        end
    end
    rows(ii)=S;
end
tableOut=struct2table(rows);
end

function S = empty_sensitivity(methods)
S=struct('caseNumber',NaN,'caseLabel',"",'discharge_Q_m3ps',NaN, ...
    'nCommonValidTransects',0,'khMedian',NaN,'khQ25',NaN,'khQ75',NaN, ...
    'khErrLow',NaN,'khErrHigh',NaN);
names={'nPairs','differenceMedian_pct','differenceQ25_pct', ...
    'differenceQ75_pct','differenceErrLow_pct','differenceErrHigh_pct'};
for ii=1:numel(methods)
    key=methods(ii).key;
    for jj=1:numel(names)
        value=NaN; if strcmp(names{jj},'nPairs'), value=0; end
        S.([key '_' names{jj}])=value;
    end
end
end
