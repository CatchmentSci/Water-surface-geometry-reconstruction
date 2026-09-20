function fig = plot_Us_vs_SWW(opts)
%% Plot surface velocity as a function of standing-wave wavelength
%
% This function calculates and plots the theoretical relationship between
% surface velocity U_s and standing-wave wavelength lambda for several
% fixed water depths and three assumed vertical velocity profiles:
% constant, linear, and power-law.
%
% The deep-water relationship U_s = sqrt(g/k) is shown for reference.
% The secondary y-axis identifies the curves corresponding to the
% different water depths.
%
% Inputs:
%   opts  Optional structure controlling figure output and appearance.
%
% Principal opts fields:
%   alpha           velocity-index (default: 0.85)
%   figureWidthIn   figure width in inches (default: 8)
%   figureHeightIn  figure height in inches (default: 3)
%   saveFigure      save figure when true (default: false)
%   outFigureFile   extensionless output filename
%   tickFontSize    source tick font size for 0.6\textwidth (default: 18)
%   labelFontSize   source label size for 0.6\textwidth (default: 20)
%   legendFontSize  source legend size for 0.6\textwidth (default: 18)
%
% Output:
%   fig   MATLAB figure handle.
%
% Example:
%   fig = plot_Us_vs_SWW();
%
% Author: Giulio Dolcetti
% Date: 16/09/2026
%
% -------------------------------------------------------------------------

if nargin < 1 || isempty(opts)
    opts = struct;
end

opts = set_plot_defaults(opts);

%% Physical parameters

g = 9.81;

%% Define wavelength and water-depth ranges

% Wavelengths considered in the figure
lambda = logspace(-1, 2, 100);

% Fixed water depths used to generate the curves
h = logspace(-1, 1, 5);

% Wavenumber
[L, H] = meshgrid(lambda, h);
K = 2 * pi ./ L;

%% Calculate surface velocity

alpha = opts.alpha;

% Constant velocity profile
U_const = sqrt( ...
    g .* H .* tanh(K .* H) ./ (K .* H));

% Linear velocity profile
m = 2 * (1 - alpha);

U_lin = sqrt( ...
    g .* H .* tanh(K .* H) ./ ...
    (K .* H - m .* tanh(K .* H)));

% Power-law velocity profile
n = (1 - alpha) / alpha;
s = sign(0.5 - n);

nu_num = s * (0.5 - n);
nu_den = -s * (0.5 + n);

U_pow = sqrt( ...
    g .* H .* ...
    besseli(nu_num, K .* H) ./ ...
    besseli(nu_den, K .* H) ./ ...
    (K .* H));

%% Define curve colours

constantColour = [0 0 0];
linearColour = [0 0 1];
powerColour = [255 127 0] / 255;

% The figure is placed at 0.6\textwidth. Its exported bounding box includes
% the right-hand depth labels, giving a measured final scale of about 0.5.
% Pre-scale visual elements so their manuscript sizes match other figures.
latexScale = 0.5;
profileLineWidth = 1.0 / latexScale;
deepLineWidth = 0.7 / latexScale;
axesLineWidth = 0.8 / latexScale;

%% Create figure

% Size the exported artwork for its full-width (5.5-inch) placement in
% the AGU manuscript, avoiding enlargement of text and strokes by LaTeX.
fig = figure('Position',[100 100 928 530]);

ax1 = axes(fig);
hold(ax1, 'on');

%% Plot velocity-profile relationships

% Store one handle for each profile family. These representative curves
% are used to construct the legend; all depth-specific curves are plotted.
for i = 1:length(h)

    p1{i}=plot(ax1, lambda, U_const(i, :), ...
        'Color', constantColour, ...
        'LineWidth', profileLineWidth);

    p2{i}=plot(ax1, lambda, U_lin(i, :), ...
        'Color', linearColour, ...
        'LineWidth', profileLineWidth);

    p3{i}=plot(ax1, lambda, U_pow(i, :), ...
        'Color', powerColour, ...
        'LineWidth', profileLineWidth);

end

% Deep-water relationship
U_deep = sqrt(g / (2 * pi) .* lambda);

p0 = plot(ax1, lambda, U_deep, ...
    'k--', ...
    'LineWidth', deepLineWidth);

%% Configure primary axes


set(ax1, ...
    'XScale', 'log', ...
    'YScale', 'log', ...
    'FontSize', opts.tickFontSize, ...
    'LineWidth', axesLineWidth, ...
    'TickDir', 'out', ...
    'Layer', 'top');

ax1.XTick = [1e-1 1e0 1e1 1e2];


xlim(ax1, [1e-1 1e2]);
ylim(ax1, [0.4 10]);

xlabel(ax1, '$\lambda$ ($\mathrm{m}$)', ...
    'Interpreter', 'latex', 'FontSize', opts.labelFontSize);

ylabel(ax1, '$U_s$ ($\mathrm{m}\,\mathrm{s}^{-1}$)', ...
    'Interpreter', 'latex', 'FontSize', opts.labelFontSize);

%% Legend

lgd = legend([p0, p1{2}, p2{2}, p3{2}], ...
    {'Deep-water limit', 'Constant profile', ...
     'Linear profile', 'Power profile'}, ...
    'Orientation', 'vertical', ...
    'Location', 'northwest', ...
    'Interpreter', 'none', ...
    'FontSize', opts.legendFontSize, ...
    'Box', 'off');
lgd.Units = 'normalized';

%% Create secondary axes for water-depth labels

% The secondary y-axis uses the same coordinates as the primary axis.
% Its tick positions correspond to the surface velocities at the largest
% wavelength for each of the five water depths.
ax2 = axes( ...
    'Position', ax1.Position, ...
    'Color', 'none', ...
    'XAxisLocation', 'top', ...
    'YAxisLocation', 'right', ...
    'XColor', 'none', ...
    'YColor', 'k', ...
    'XScale', 'log', ...
    'YScale', 'log', ...
    'FontSize', opts.tickFontSize, ...
    'LineWidth', axesLineWidth, ...
    'TickDir', 'out', ...
    'Layer', 'top');

linkaxes([ax1 ax2], 'xy');

% Use the constant-profile velocity at the largest wavelength as the
% representative velocity associated with each water depth.
yTicks = U_const(:, end);

ax2.YTick = yTicks;
ax2.YTickLabel = compose('$h=%.1f~\\mathrm{m}$', h);
ax2.FontName = 'Times New Roman';
ax2.TickLabelInterpreter = 'latex';

% Keep the two axes exactly coincident.
ax1.Position = [0.12 0.2 0.45 0.7];
ax2.Position = ax1.Position;

ylim(ax1, [0.4 10]);

% Ensure the secondary axis does not alter the displayed limits.
xlim(ax2, ax1.XLim);
ylim(ax2, ax1.YLim);


%% Figure appearance

set(fig, 'Color', 'w');


%% Save figure

if opts.saveFigure && ~isempty(opts.outFigureFile)

    [outDir, ~, ext] = fileparts(opts.outFigureFile);

    if ~isempty(outDir) && ~isfolder(outDir)
        mkdir(outDir);
    end

    if isempty(ext)

        pngFile = [opts.outFigureFile '.png'];
        pdfFile = [opts.outFigureFile '.pdf'];

        exportgraphics(fig, pngFile, ...
            'Resolution', 600);

        exportgraphics(fig, pdfFile, ...
            'ContentType', 'vector', ...
            'BackgroundColor', 'white');

    else

        exportgraphics(fig, opts.outFigureFile, ...
            'Resolution', 600);

    end

end

end


%% ========================================================================
% Local helpers
%% ========================================================================

function opts = set_plot_defaults(opts)

opts = default_field(opts, 'alpha', 0.85);
opts = default_field(opts, 'figureWidthIn', 8);
opts = default_field(opts, 'figureHeightIn', 3);
opts = default_field(opts, 'saveFigure', false);
opts = default_field(opts, 'outFigureFile', '');
opts = default_field(opts, 'tickFontSize', 9 / 0.5);
opts = default_field(opts, 'labelFontSize', 10 / 0.5);
opts = default_field(opts, 'legendFontSize', 9 / 0.5);

end


function s = default_field(s, f, v)

if ~isfield(s, f) || isempty(s.(f))
    s.(f) = v;
end

end
