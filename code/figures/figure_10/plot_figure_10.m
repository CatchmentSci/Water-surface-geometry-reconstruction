function outputs = plot_figure_10(opts)
%PLOT_FIGURE_10 Compare predicted depth based on empirical relationship with observed one.
%
% This function reproduces Figure 10, comparing observed depth with
% prediction based on water-surface wavelength and amplitude.
%
% Required input files:
%   per_transect_initial_accepted.csv
%   Dart_video_statistics.xlsx
%   videos/derived/profiles/*.csv
%
% Required function:
%   Fr_calc.m
%
% Outputs:
%   outputs - structure containing the calculated quantities and
%             paths to the generated figures.

%% Defaults

if ~isfield(opts, 'dataRoot')
    opts.dataRoot = "";
end

if ~isfield(opts, 'saveFigure')
    opts.saveFigure = false;
end

if ~isfield(opts, 'outFigureFile')
    opts.outFigureFile = '';
end


%% Parameters

g = 9.81;
alpha = 0.85;
ks = 0.05;
lambdaMinThreshold = 1;
lambdaMaxThreshold = 7;

org = [255,127,0] / 255;
cmap = local_coolwarm(256);

% The AGU manuscript displays this 4-inch-wide figure at 0.6 of its
% 5.5-inch text width. Scale the artwork so its final sizes match the
% full-width figures (9/10-pt text and 30-pt^2 markers).
latexScale = (0.6 * 5.5) / 4.0;
axesFontSize = 9 / latexScale;
labelFontSize = 10 / latexScale;
colorbarFontSize = 9 / latexScale;
MS = sqrt(30) / latexScale;
axesLineWidth = 0.8 / latexScale;
dataLineWidth = 0.7 / latexScale;
errorColour = [0.62 0.62 0.62];
markerEdgeColour = [0.15 0.15 0.15];

%% Input files

inputDataFolder = char(opts.dataRoot);

velocityTableFile = fullfile( ...
    inputDataFolder, 'per_transect_initial_accepted.csv');

if ~isfile(velocityTableFile)
    matches = dir(fullfile( ...
        inputDataFolder, '**', 'per_transect_initial_accepted.csv'));
    if numel(matches) ~= 1
        error('Expected one per_transect_initial_accepted.csv below %s.', ...
            inputDataFolder);
    end
    velocityTableFile = fullfile(matches(1).folder, matches(1).name);
    inputDataFolder = matches(1).folder;
end

dataFolderProfiles = fullfile( ...
    inputDataFolder, 'profiles');

if ~isfolder(dataFolderProfiles)
    profileMatches = dir(fullfile(char(opts.dataRoot), '**', ...
        'profiles', '*.csv'));
    profileFolders = unique(string({profileMatches.folder}));
    if numel(profileFolders) ~= 1
        error('Expected one derived profiles folder below %s.', ...
            opts.dataRoot);
    end
    dataFolderProfiles = char(profileFolders(1));
end

dsv = readtable(velocityTableFile);

stats = readtable(fullfile( ...
    inputDataFolder, 'Dart_video_statistics.xlsx'));

%% Sort observations by discharge

[~, sortQindx] = sort(stats.Discharge_m3s, 'descend');
stats = stats(sortQindx, :);

dates = datetime(stats.VideoDatetime);
Q = stats.Discharge_m3s;

Qmin = min(Q);
Qmax = max(Q);

idxColor = round(1 + ...
    (Q - Qmin) ./ (Qmax - Qmin) * 255);

QColor = cmap(idxColor, :);

%% Process observations

nCases = length(dates);

csVel = cell(nCases, 1);
csDepth = cell(nCases, 1);
csFr = cell(nCases, 1);
csLambda = cell(nCases, 1);
csAmpl = cell(nCases, 1);
csKd = cell(nCases, 1);

velMean = nan(nCases, 1);
velMin = nan(nCases, 1);
velMax = nan(nCases, 1);

depthMean = nan(nCases, 1);
depthMin = nan(nCases, 1);
depthMax = nan(nCases, 1);

FrMean = nan(nCases, 1);
FrMin = nan(nCases, 1);
FrMax = nan(nCases, 1);

lambdaMean = nan(nCases, 1);
lambdaMin = nan(nCases, 1);
lambdaMax = nan(nCases, 1);

amplMean = nan(nCases, 1);
amplMin = nan(nCases, 1);
amplMax = nan(nCases, 1);

kdMean = nan(nCases, 1);
kdMin = nan(nCases, 1);
kdMax = nan(nCases, 1);

for id = 1:nCases

    % Surface velocity
    velIndx = find(dsv.xCase == "R" + id);

    if ~isempty(velIndx)
        vel = dsv.U_initial(velIndx);
        csVel{id} = vel;
    else
        csVel{id} = nan;
    end

    velMean(id) = nanmedian(csVel{id});
    velMin(id) = prctile(csVel{id}, 25);
    velMax(id) = prctile(csVel{id}, 75);

    % Profile data
    filename = char(string( ...
        dates(id), 'yyyyMMdd_HHmms'));

    ds = readtable(fullfile( ...
        dataFolderProfiles, ...
        ['devon_dart', filename, ...
         '_selected_map_profiles_for_real_batch_summary.csv']));

    % Depth
    idx = find(~strcmp(ds.crossSectionDepthStatus, 'ok'));

    depth = ds.crossSectionDepth_m;
    depth(idx) = NaN;

    csDepth{id} = depth;

    depthMean(id) = nanmedian(csDepth{id});
    depthMin(id) = prctile(csDepth{id}, 25);
    depthMax(id) = prctile(csDepth{id}, 75);

    % Froude number
    csFr{id} = vel ./ sqrt(g * depth);

    FrMean(id) = nanmedian(csFr{id});
    FrMin(id) = prctile(csFr{id}, 25);
    FrMax(id) = prctile(csFr{id}, 75);

    % Wavelength
    lambda = accepted_autocorr_wavelengths( ...
        ds.autocorrWavelength_m, ds.autocorrPeakR, ...
        lambdaMinThreshold, lambdaMaxThreshold);

    csLambda{id} = lambda;

    lambdaMean(id) = nanmedian(csLambda{id});
    lambdaMin(id) = prctile(csLambda{id}, 25);
    lambdaMax(id) = prctile(csLambda{id}, 75);

    % Wave amplitude
    ampl = ds.amplitude98_2_m;
    flagRetain = ds.amplitudeRetainedForSummary;

    ampl(flagRetain == 0) = NaN;
    csAmpl{id} = ampl;

    amplMean(id) = nanmedian(csAmpl{id});
    amplMin(id) = prctile(csAmpl{id}, 25);
    amplMax(id) = prctile(csAmpl{id}, 75);

    % kh
    csKd{id} = 2 * pi ./ lambda .* depth;

    kdMean(id) = nanmedian(csKd{id});
    kdMin(id) = prctile(csKd{id}, 25);
    kdMax(id) = prctile(csKd{id}, 75);

end

%% Fit amplitude-depth relationship using all data

idx = 1:nCases;
[pfitAllData, ~] = polyfit( ...
    log(kdMean(idx)), ...
    log(amplMean(idx)), ...
    1);


%% cross-validation approach using 3 points at a time

LocalErr = []; %cell(nCases,1);
LocalErrRelative = [];

dEstCv = cell(nCases,1);
dEstCvMean = nan(nCases,1);
dEstCvMin = nan(nCases,1);
dEstCvMax = nan(nCases,1);

n_points = 3; % number of points used for a single calibration
C = nchoosek(idx,n_points); % create combinations of 3 points
n = size(C,1);

c1 = nan(n,1);
c2 = nan(n,1);

for i=1:size(C,1)
    idx_train = C(i,:);                 % training set
    idx_test = setdiff(idx,idx_train);  % testing set

    % fit for a 3 cases combination
    pfit = polyfit(log(kdMean(idx_train)),log(amplMean(idx_train)),1);

    c1(i) = exp(pfit(2));   % coeff 1
    c2(i) = pfit(1);        % coeff 2

    for ii=idx_test
    % invert to estimate depth of testing set
        pfitInv = (csAmpl{ii}/c1(i)) .^(1/c2(i)) .* csLambda{ii} / (2*pi);
        pfitInv(pfitInv > 10) = NaN; % rejects unreasonably high depth estimations > 10 m
        dEstCv{ii} = [dEstCv{ii} ; pfitInv];
        LocalErr = [LocalErr; pfitInv - csDepth{ii}];
        LocalErrRelative = [LocalErrRelative; (pfitInv - csDepth{ii})./csDepth{ii}];
    end

end

% calculate statitics of reconstructions
for i = 1:nCases
dEstCvMean(i) = nanmedian(dEstCv{i},1);
dEstCvMin(i) = prctile(dEstCv{i},25);
dEstCvMax(i) = prctile(dEstCv{i},75);
end

% statistics of fitted parameters
c1Mean = nanmedian(c1);
c1Min = prctile(c1,25);
c1Max = prctile(c1,75);

c2Mean = nanmedian(c2);
c2Min = prctile(c2,25);
c2Max = prctile(c2,75);

% coefficient of determination
SSRes = sum((dEstCvMean - depthMean).^2);
SSTot = sum((dEstCvMean - mean(dEstCvMean)).^2);
R2 = 1 - SSRes/SSTot;

% root-mean-squared error
RMS = sqrt(mean((dEstCvMean - depthMean).^2));
RMSRelative =   sqrt(mean((dEstCvMean - depthMean).^2./depthMean.^2));

fprintf('\nr2 depth estimation (median): %.2f',R2)
fprintf('\nrms depth estimation (median): %.2f',RMS)
fprintf('\nrms relative depth estimation (median): %.2f',RMSRelative)

LocalIR = mean(dEstCvMax-dEstCvMin)/2;
LocalIRRelative = mean((dEstCvMax-dEstCvMin)./dEstCvMean)/2;

fprintf('\n average interquartile range for local depth reconstruction: %.2f',LocalIR)
fprintf('\n average relative interquartile range for local depth reconstruction: %.2f',LocalIRRelative)
fprintf('\n average measurement uncertainty: %.2f',mean(depthMax-depthMin)/2)

LocalErrMedian = nanmedian(abs(LocalErr));
LocalErrRelativeMedian = nanmedian(abs(LocalErrRelative));
fprintf('\n median absolute error of local depth reconstruction: %.2f',LocalErrMedian)
fprintf('\n median relative error of local depth reconstruction: %.2f',LocalErrRelativeMedian)


%% Figure 10: observed vs estimated depth

fg = figure( ...
    'Units', 'inches', ...
    'Position', [1 1 4.0 3.15], ...
    'PaperPositionMode', 'auto', ...
    'Color', 'w');

ax = axes(fg);
hold(ax, 'on');

plot(ax, [1 3], [1 3], '--', ...
    'Color', [0.20 0.20 0.20], ...
    'LineWidth', axesLineWidth, ...
    'HandleVisibility', 'off');

for id = 1:nCases

    errorbar(ax, ...
        depthMean(id), dEstCvMean(id), ...
        dEstCvMean(id) - dEstCvMin(id), ...
        dEstCvMax(id) - dEstCvMean(id), ...
        depthMean(id) - depthMin(id), ...
        depthMax(id) - depthMean(id), ...
        'o', ...
        'Color', errorColour, ...
        'LineStyle', 'none', ...
        'MarkerSize', MS, ...
        'MarkerEdgeColor', markerEdgeColour, ...
        'MarkerFaceColor', QColor(id, :), ...
        'LineWidth', dataLineWidth);

end

xlabel(ax, '$h$ ($\mathrm{m}$)', ...
    'Interpreter', 'latex', 'FontSize', labelFontSize);
ylabel(ax, '$h_{\mathrm{est}}$ ($\mathrm{m}$)', ...
    'Interpreter', 'latex', 'FontSize', labelFontSize);

colormap(ax, cmap);
clim(ax, [Qmin Qmax]);

cb = colorbar(ax, 'eastoutside');
cb.Label.Interpreter = 'latex';
cb.Label.String = '$Q\;(\mathrm{m}^{3}\,\mathrm{s}^{-1})$';
cb.Label.FontSize = labelFontSize;
cb.FontSize = colorbarFontSize;
cb.TickLabelInterpreter = 'tex';
cb.TickDirection = 'out';

xlim(ax, [1 3]);
ylim(ax, [1 3]);
xticks(ax, 1:1:3);
yticks(ax, 1:1:3);
axis(ax, 'square');
box(ax, 'on');
set(ax, 'FontSize', axesFontSize, 'LineWidth', axesLineWidth, ...
    'TickDir', 'out', 'TickLabelInterpreter', 'tex', 'Layer', 'top');

if opts.saveFigure

    exportgraphics(fg, ...
        [opts.outFigureFile '.png'], ...
        'ContentType', 'vector', ...
        'Resolution', 600);

    exportgraphics(fg, ...
        [opts.outFigureFile '.pdf'], ...
        'ContentType', 'vector', ...
        'Resolution', 600);

end


%% Outputs

outputs = struct();

outputs.depth = csDepth;
outputs.velocity = csVel;
outputs.froude = csFr;
outputs.wavelength = csLambda;
outputs.amplitude = csAmpl;
outputs.kh = csKd;

outputs.depthEstimated = dEstCv;

outputs.depthMean = depthMean;
outputs.depthEstimatedMean = dEstCvMean;
outputs.velocityMean = velMean;
outputs.froudeMean = FrMean;
outputs.wavelengthMean = lambdaMean;
outputs.amplitudeMean = amplMean;
outputs.khMean = kdMean;

outputs.amplitudeFit = pfitAllData;
outputs.amplitudeFitCv = [c1'; c2'];

outputs.R2 = R2;
outputs.rmsError = RMS;

outputs.LocalIR = LocalIR;
outputs.LocalIRRelative = LocalIRRelative;

outputs.LocalErrMedian = LocalErrMedian;
outputs.LocalErrRelativeMedian = LocalErrRelativeMedian;

outputs.figure10 = opts.outFigureFile;


end


%% Local functions

function add_panel_label(ax, label, fontSize)

text(ax, -0.3, 1.0, label, ...
    'Units', 'normalized', ...
    'FontSize', fontSize, ...
    'Interpreter', 'latex');

end

function map = local_coolwarm(m)
if nargin < 1 || isempty(m)
    m = 256;
end
if m <= 0
    map = zeros(0, 3);
    return
end
anchors = [59 76 192; 84 112 222; 129 164 251; 180 205 251; ...
    221 221 221; 241 184 156; 229 112 88; 203 62 56; 180 4 38] ./ 255;
map = interp1(linspace(0, 1, size(anchors, 1)), anchors, ...
    linspace(0, 1, m), "linear");
map = max(0, min(1, map));
end


function lambda = accepted_autocorr_wavelengths( ...
        lambda, peakR, lambdaMin, lambdaMax)
% Apply the canonical field-case wavelength quality and outlier filters.
lambda = double(lambda(:));
peakR = double(peakR(:));
lambda(~isfinite(lambda) | ~isfinite(peakR) | peakR < 0.10 | ...
    lambda < lambdaMin | lambda > lambdaMax) = NaN;
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
