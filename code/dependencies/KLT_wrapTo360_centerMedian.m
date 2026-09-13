function [ang_out, shift] = KLT_wrapTo360_centerMedian(ang)
% ============================================================
% Shift angles so that their circular mean is at 180 degrees
% and wrap output to [0, 360)
%
% INPUT:
%   ang      - input angles (degrees, any range)
%
% OUTPUT:
%   ang_out  - angles in [0, 360)
%   shift    - scalar shift applied (degrees)
% ============================================================

% store original shape
sz = size(ang);

% reshape to column for processing
ang = ang(:);

% convert to radians
ang_rad = deg2rad(ang);

% circular mean
mean_ang = atan2(mean(sin(ang_rad)), mean(cos(ang_rad)));

% convert to degrees
mean_ang_deg = rad2deg(mean_ang);

% compute required shift to move mean -> 180
shift = 180 - mean_ang_deg;

% apply shift
ang_shifted = ang + shift;

% wrap to [0, 360)
ang_out = mod(ang_shifted, 360);

% restore original shape
ang_out = reshape(ang_out, sz);

end