function map = figure_01_colormap(m)
%FIGURE_01_COLORMAP Cool-to-warm diverging colour map for Figure 1.
%   This small repository-owned implementation interpolates between chosen
%   cool, neutral and warm RGB anchors and avoids an external colormap file.

arguments
    m (1,1) double {mustBeInteger, mustBeNonnegative} = 256
end

if m == 0
    map = zeros(0, 3);
    return
end

anchors = [59 76 192; 221 221 221; 180 4 38] ./ 255;
anchorPositions = [0 0.5 1];
map = interp1(anchorPositions, anchors, linspace(0, 1, m), 'pchip');
map = max(0, min(1, map));
end
