function fig = plot_kh_vs_F(opts)
%% Plot Froude number as a function of kh for different velocity profiles
%
% This script calculates and plots the relationship between the
% dimensionless depth kh and the Froude number (Fr) for three
% vertical velocity profiles:
%   1. constant velocity
%   2. linear velocity profile
%   3. power-law velocity profile
%
% The second panel shows the corresponding
% relationship scaled by alpha, where alpha is the velocity-index.
%
% The dashed line represents the deep-water approximation:
%   Fr = (kh)^(-1/2)
%
% Inputs:
%   opts  optional structure controlling figure output and appearance.
%
% Principal opts fields:
%   alpha          velocity-index (default: 0.85)
%   figureWidthIn  figure width in inches (default: 8)
%   figureHeightIn figure height in inches (default: 3)
%   saveFigure     save figure when true (default: false)
%   outFigureFile  extensionless output filename
%
% Output:
%   fig   MATLAB figure handle.
%
% Example:
%   fig = plot_kh_vs_F();
%
% Author: Giulio Dolcetti
% Date: 16/09/2026
%
% -------------------------------------------------------------------------

if nargin < 1 || isempty(opts)
    opts = struct;
end
opts = set_plot_defaults(opts);

%% Define dimensionless depth values

kh = logspace(-1,2,100);

%% Velocity index
alpha = opts.alpha;

%% Constant velocity profile
%  Fr^2 = tanh(kh)/(kh).

F2_const = tanh(kh)./kh;

%% Linear velocity profile
% 

m = 2*(1-alpha);

F2_lin = tanh(kh)./(kh - m*tanh(kh));

%% Power-law velocity profile
%

n = 1/alpha - 1;
s = sign(0.5 - n);

F2_pow = besseli(s*(0.5-n),kh) ./ ...
         besseli(-s*(0.5+n),kh) ./ kh;

%% Plot

org = [255,127,0]/255;

fig = figure('Color','w', ...
        'Units','inches', ...
        'Position',[0.5 0.5 opts.figureWidthIn opts.figureHeightIn], ...
        'PaperUnits','inches', ...
        'PaperPosition',[0 0 opts.figureWidthIn opts.figureHeightIn], ...
        'PaperSize',[opts.figureWidthIn opts.figureHeightIn], ...
        'InvertHardcopy','off');

    tl = tiledlayout(fig,1,2, ...
        'TileSpacing','compact', ...
        'Padding','compact');

%% Panel (a): Froude number

ax1 = nexttile;
hold on

plot(kh,sqrt(F2_const),'k','LineWidth',1.5)
plot(kh,sqrt(F2_lin),'b','LineWidth',1.5)
plot(kh,sqrt(F2_pow),'color',org,'LineWidth',1.5)
plot(kh,1./kh.^(1/2),'k--','LineWidth',1.5)

set(ax1,'XScale','log')

xlabel('$kh$','Interpreter','latex');
ylabel('$\rm{Fr}$','Interpreter','latex');

set(ax1,'XScale','log', ...
        'XTick',[1e-1 1e0 1e1 1e2], ...
        'FontSize',opts.tickFontSize);

xlabel(ax1,'$kh$','Interpreter','latex');
ylabel(ax1,'$Fr$','Interpreter','latex');
ylim(ax1,[0 1.5]);

add_panel_label(ax1,'(a)',opts);


%% Panel (b): alpha-scaled Froude number

ax2 = nexttile(tl,2);
hold(ax2,'on');

p1 = plot(ax2,kh,sqrt(F2_const),'k','LineWidth',1.5);
p2 = plot(ax2,kh,alpha*sqrt(F2_lin),'b','LineWidth',1.5);
p3 = plot(ax2,kh,alpha*sqrt(F2_pow),'Color',org,'LineWidth',1.5);
p4 = plot(ax2,kh,1./sqrt(kh),'k--','LineWidth',1.5);

set(ax2,'XScale','log', ...
    'XTick',[1e-1 1e0 1e1 1e2], ...
    'FontSize',opts.tickFontSize);

xlabel(ax2,'$kh$','Interpreter','latex');
ylabel(ax2,'$\alpha Fr$','Interpreter','latex');
ylim(ax2,[0 1.5]);

add_panel_label(ax2,'(b)',opts);


%% Legend

lgd = legend(ax2,[p1 p2 p3 p4], ...
    {'constant','linear','power-funct.','deep water'}, ...
    'Orientation','vertical', ...
    'Box','off');

lgd.Layout.Tile = 'east';
lgd.FontSize = opts.legendFontSize;

%% Save figure

if opts.saveFigure && ~isempty(opts.outFigureFile)
    [outDir,~,ext] = fileparts(opts.outFigureFile);

    if ~isempty(outDir) && ~isfolder(outDir)
        mkdir(outDir);
    end

    if isempty(ext)
        pngFile = [opts.outFigureFile '.png'];
        pdfFile = [opts.outFigureFile '.pdf'];

        exportgraphics(fig,pngFile,'Resolution',600);
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
    opts = default_field(opts,'tickFontSize',16);
    opts = default_field(opts,'legendFontSize',16);
    opts = default_field(opts,'panelLabelFontSize',16);
end


function s = default_field(s,f,v)

    if ~isfield(s,f) || isempty(s.(f))
        s.(f) = v;
    end
end


function add_panel_label(ax,txt,opts)

    text(ax,0.9,0.98,txt, ...
        'Units','normalized', ...
        'FontSize',opts.panelLabelFontSize, ...
        'VerticalAlignment','top', ...
        'HorizontalAlignment','left', ...
        'Color','k');
end