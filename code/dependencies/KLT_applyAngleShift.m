function ang_out = KLT_applyAngleShift(ang, shift)
% ============================================================
% Apply a given angular shift and wrap to [0, 360)
%
% INPUT:
%   ang      - input angles in degrees (any shape)
%   shift    - scalar shift in degrees
%
% OUTPUT:
%   ang_out  - shifted angles in [0, 360)
%
% NOTES:
%   - Preserves NaN values
% ============================================================

ang_out = mod(ang + shift, 360);

% preserve NaNs from input
ang_out(~isfinite(ang)) = NaN;

end