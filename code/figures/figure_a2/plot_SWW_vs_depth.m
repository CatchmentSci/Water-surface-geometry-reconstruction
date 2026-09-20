function fig = plot_SWW_vs_depth(opts)
%% Plot standing-wave wavelength as a function of water depth
%
% This function calculates and plots the theoretical relationship between
% standing-wave wavelength lambda and water depth h for several
% surface flow velocities and three assumed vertical velocity profiles:
% constant, linear, and power-law.
%
% Lines of constant Froude number are shown for reference.
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
%   tickFontSize    axis tick font size (default: 13)
%   labelFontSize   axis-label font size (default: 14)
%   legendFontSize  legend font size (default: 12)
%
% Output:
%   fig   MATLAB figure handle.
%
% Example:
%   fig = plot_SWW_vs_depth();
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
alpha = opts.alpha;

%% Calculate wavelength as a function of depth

h = logspace(-1,1,100);
u = logspace(0,1,10);

[H,UH] = meshgrid(h,u);

l_const = zeros(size(H));
l_lin   = zeros(size(H));
l_pow   = zeros(size(H));

for i = 1:length(u)
    l_const(i,:) = expected_wavelength( ...
        H(i,:),UH(i,:),1,'constant');

    l_lin(i,:) = expected_wavelength( ...
        H(i,:),UH(i,:),alpha,'linear');

    l_pow(i,:) = expected_wavelength( ...
        H(i,:),UH(i,:),alpha,'power');
end

% Logarithmic axes cannot display nonpositive failed-root values. Converting
% only those values to NaN preserves the calculations and visible curves
% while avoiding MATLAB's "Negative data ignored" graphics warning.
l_const(l_const <= 0) = NaN;
l_lin(l_lin <= 0) = NaN;
l_pow(l_pow <= 0) = NaN;

%% Reconstruct depth as a function of wavelength

l = logspace(-1,1,100);
u = logspace(0,1,10);

[L,UL] = meshgrid(l,u);

h_const = zeros(size(L));
h_lin   = zeros(size(L));
h_pow   = zeros(size(L));

for i = 1:length(u)

    dmax = max(UL(i,:).^2 / (g*0.8^2));

    h_const(i,:) = reconstr_depth( ...
        L(i,:),UL(i,:),1,[1e-01 dmax],'constant');

    h_lin(i,:) = reconstr_depth( ...
        L(i,:),UL(i,:),alpha,[1e-01 dmax],'linear');

    h_pow(i,:) = reconstr_depth( ...
        L(i,:),UL(i,:),alpha,[1e-01 dmax],'power');

end


h_const(h_const <= 0) = NaN;
h_lin(h_lin <= 0) = NaN;
h_pow(h_pow <= 0) = NaN;

%% Froude-number reference lines

Fr = [0.3 0.4 0.55 0.75 0.90 0.97 0.99];

kh = expected_kh(Fr,1,'constant');

%% Colours

constantColour = [0 0 0];
linearColour = [0 0 1];
powerColour = [255 127 0] / 255;
froudeColour = [0.55 0.55 0.55];

%% Create figure

fig = figure('Position',[100 100 700 450]);

%% Plot Froude-number reference lines

for i = 1:length(Fr)
    hold on
    plot(2*pi*h./kh(i),h, ...
        '--', ...
        'Color',froudeColour, ...
        'LineWidth',0.7);
end

%% Plot theoretical relationships

for i = 1:length(u)

    % Forward calculation: wavelength as a function of depth
    hold on
    p1{i} = plot(l_const(i,:),H(i,:), ...
        'Color',constantColour, ...
        'LineWidth',1.2);

    hold on
    p2{i} = plot(l_lin(i,:),H(i,:), ...
        'Color',linearColour, ...
        'LineWidth',1.2);

    hold on
    p3{i} = plot(l_pow(i,:),H(i,:), ...
        'Color',powerColour, ...
        'LineWidth',1.2);

    % Inverse calculation: reconstructed depth as a function of wavelength
    hold on
    plot(l,h_const(i,:), ...
        'Color',constantColour, ...
        'LineWidth',1.2);

    hold on
    plot(l,h_lin(i,:), ...
        'Color',linearColour, ...
        'LineWidth',1.2);

    hold on
    plot(l,h_pow(i,:), ...
        'Color',powerColour, ...
        'LineWidth',1.2);

end

%% Axes

ax1 = gca;

set(ax1, ...
    'XScale','log', ...
    'YScale','log', ...
    'FontSize',opts.tickFontSize, ...
    'LineWidth',0.8, ...
    'TickDir','out', ...
    'TickLabelInterpreter','tex', ...
    'Layer','top');
box(ax1,'on');

lgd = legend([p1{3},p2{3},p3{3}], ...
    {'Constant profile','Linear profile','Power profile'}, ...
    'Orientation','vertical', ...
    'Location','eastoutside', ...
    'FontSize',opts.legendFontSize, ...
    'Interpreter','none', ...
    'Box','off');

lgd.Units = 'normalized';

xlabel(ax1,'$\lambda$ ($\mathrm{m}$)', ...
    'Interpreter','latex','FontSize',opts.labelFontSize);
ylabel(ax1,'$h$ ($\mathrm{m}$)', ...
    'Interpreter','latex','FontSize',opts.labelFontSize);

%% Froude-number labels

yline_array = logspace(-0.9,0.7,length(Fr));
yline_array = fliplr(yline_array);

for i = 1:length(Fr)

    yline = yline_array(i);
    xline = 2*pi*yline./kh(i);

    text(ax1,xline,yline,sprintf('$\\mathrm{Fr}=%.2f$',Fr(i)), ...
        'Interpreter','latex', ...
        'Rotation',25, ...
        'FontSize',opts.annotationFontSize, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'BackgroundColor','w', ...
        'Margin',0.1);

end

%% Secondary x-axis showing surface velocity

ax2 = axes('Position',ax1.Position, ...
    'Color','none', ...
    'XLim',ax1.XLim, ...
    'YLim',ax1.YLim, ...
    'XAxisLocation','top', ...
    'YAxisLocation','right', ...
    'YColor','none', ...
    'YScale','log', ...
    'XScale','log', ...
    'FontSize',opts.tickFontSize, ...
    'LineWidth',0.8, ...
    'TickDir','out', ...
    'TickLabelInterpreter','latex');

linkaxes([ax1 ax2],'xy');

I = find(l_const(:,end)>1,1,'last');

ax2.XTick = l_const(1:I,end);
ax2.XTickLabel = compose( ...
    '$U_s=%.1f~\\mathrm{m}\\,\\mathrm{s}^{-1}$',u(1:I));

% Preserve original axes dimensions and position.
ax1.Position = [0.12 0.15 0.55 0.62];
ax2.Position = ax1.Position;

xlim(ax1,[0.5 10]);

set(fig,'Color','w');

%% Save figure

if opts.saveFigure && ~isempty(opts.outFigureFile)

    [outDir,~,ext] = fileparts(opts.outFigureFile);

    if ~isempty(outDir) && ~isfolder(outDir)
        mkdir(outDir);
    end

    if isempty(ext)

        pngFile = [opts.outFigureFile '.png'];
        pdfFile = [opts.outFigureFile '.pdf'];

        exportgraphics(fig,pngFile, ...
            'Resolution',600);

        exportgraphics(fig,pdfFile, ...
            'ContentType','vector', ...
            'BackgroundColor','white');

    else

        exportgraphics(fig,opts.outFigureFile, ...
            'Resolution',600);

    end

end

end


%% ========================================================================
% Local helpers
%% ========================================================================

function opts = set_plot_defaults(opts)

opts = default_field(opts,'alpha',0.85);
opts = default_field(opts,'figureWidthIn',8);
opts = default_field(opts,'figureHeightIn',3);
opts = default_field(opts,'saveFigure',false);
opts = default_field(opts,'outFigureFile','');
opts = default_field(opts,'tickFontSize',13);
opts = default_field(opts,'labelFontSize',14);
opts = default_field(opts,'legendFontSize',12);
opts = default_field(opts,'annotationFontSize',10.5);

end


function s = default_field(s,f,v)

if ~isfield(s,f) || isempty(s.(f))
    s.(f) = v;
end

end
