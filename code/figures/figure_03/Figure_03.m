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
% Repository:
% https://github.com/CatchmentSci/Water-surface-geometry-reconstruction
%
% Associated publication:
% Perks & Dolcetti, Water-surface geometry reconstruction for non-contact 
% river flow monitoring through monocular imagery and inverse modeling
%
% Author: Giulio Dolcetti
% Date: 16/09/2026
%
% -------------------------------------------------------------------------

%% Define dimensionless depth values

kh = logspace(-1,2,100);

%% Velocity index
alpha = 0.85;

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

fg = figure('Position',[100 100 800 300]);
t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

%% Panel (a): Froude number

ax1 = nexttile;
hold on

plot(kh,sqrt(F2_const),'k','LineWidth',1.5)
plot(kh,sqrt(F2_lin),'b','LineWidth',1.5)
plot(kh,sqrt(F2_pow),'color',org,'LineWidth',1.5)

% Deep-water approximation
plot(kh,1./kh.^(1/2),'k--','LineWidth',1.5)

set(ax1,'XScale','log')

xlabel('$kh$','Interpreter','latex');
ylabel('$\rm{Fr}$','Interpreter','latex');

set(ax1,'XTick',[1e-1 1e0 1e1 1e2]);
set(ax1,'FontSize',16);

text(0.9,0.98,'(a)','Units','normalized','FontSize',16,...
    'VerticalAlignment','top')

ylim([0 1.5]);

%% Panel (b): alpha-scaled Froude number

ax2 = nexttile;
hold on

p1 = plot(kh,sqrt(F2_const),'k','LineWidth',1.5);
p2 = plot(kh,alpha*sqrt(F2_lin),'b','LineWidth',1.5);
p3 = plot(kh,alpha*sqrt(F2_pow),'color',org,'LineWidth',1.5);

% Deep-water approximation
p4 = plot(kh,1./kh.^(1/2),'k--','LineWidth',1.5);

set(ax2,'XScale','log')

xlabel('$kh$','Interpreter','latex');
ylabel('$\alpha \rm{Fr}$','Interpreter','latex');

set(ax2,'XTick',[1e-1 1e0 1e1 1e2]);
set(ax2,'FontSize',16);
set(ax2,'Color','w');

ylim([0 1.5]);

text(0.9,0.98,'(b)','Units','normalized','FontSize',16,...
    'VerticalAlignment','top')

% Shared legend placed outside the tiled layout
lgd = legend([p1 p2 p3 p4],...
    {'constant','linear','power-funct.','deep water'},...
    'Orientation','vertical');

lgd.Layout.Tile = 'east';

set(gcf,'Color','w');

%% Export figure

saveas(fg,'Figure_03.jpg');

exportgraphics(fg,'Figure_03.pdf',...
    'ContentType','vector',...
    'Resolution',300);