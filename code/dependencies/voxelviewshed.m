function vis = voxelviewshed(X,Y,Z,camxyz)
% Calculate a viewshed over a DEM (same logic, cached + griddedInterpolant)
%
% USAGE:
%    vis = voxelviewshed(X,Y,Z,camxyz)
%
% INPUTS:
%    X,Y,Z   : regular-grid DEM
%    camxyz  : [x y z] viewpoint
%
% OUTPUT:
%    vis     : logical visibility matrix, same size as Z

if nargin==0
    warning('no input arguments. Will use test dataset.... ')
    [X,Y] = meshgrid(-299:300);
    Z = abs(ifft2(ifftshift((hypot(Y,X)+1e-5).^(-2.1).*exp(rand(size(X))*2i*pi))));
    Z = (Z-Z(300,300))/std(Z(:));
    camxyz = [0.002,0.002,0.5];

    vis = voxelviewshed(X,Y,Z,camxyz);

    clf;
    surf(X,Y,Z,Z.*vis-~vis*2,'EdgeColor','none','FaceColor','interp');
    hold on; colormap jet; camlight;
    plot3(camxyz(1),camxyz(2),camxyz(3),'k.','markersize',40)
    plot3([0 camxyz(1)],[0 camxyz(2)],[0 camxyz(3)],'k-','linewidth',3)
    clear vis
    title('monotone=hidden from view | BlackDot=camera position')
    return
end

% ---------------------------------------------------------------------
% Options
% ---------------------------------------------------------------------
useSingle = true;   % set false to keep double precision working arrays

% ---------------------------------------------------------------------
% Persistent cache for geometry-dependent terms
% ---------------------------------------------------------------------
persistent cache

if ~isequal(size(X), size(Y), size(Z))
    error('X, Y, and Z must have the same size.')
end
if numel(camxyz) ~= 3
    error('camxyz must be a 3-element vector [x y z].')
end

sz = size(Z);

dx0 = abs(X(2,2)-X(1,1));
dy0 = abs(Y(2,2)-Y(1,1));

if useSingle
    Xw   = single(X);
    Yw   = single(Y);
    Zw   = single(Z);
    camw = single(camxyz(:).');
    dx   = single(dx0);
    dy   = single(dy0);
else
    Xw   = X;
    Yw   = Y;
    Zw   = Z;
    camw = double(camxyz(:).');
    dx   = dx0;
    dy   = dy0;
end

% ---------------------------------------------------------------------
% Rebuild cache only if geometry changed
% ---------------------------------------------------------------------
needRebuild = true;

if ~isempty(cache)
    sameSize  = isequal(cache.sz, sz);
    sameX     = isequaln(cache.X, Xw);
    sameY     = isequaln(cache.Y, Yw);
    sameCamXY = isequaln(cache.camxy, camw(1:2));
    sameDxDy  = isequaln(cache.dx, dx) && isequaln(cache.dy, dy);

    needRebuild = ~(sameSize && sameX && sameY && sameCamXY && sameDxDy);
end

if needRebuild
    Xn = (Xw(:) - camw(1)) ./ dx;
    Yn = (Yw(:) - camw(2)) ./ dy;

    rxy  = hypot(Xn, Yn);
    xang = (atan2(Yn, Xn) + pi) / (2*pi);

    rbin = round(rxy);
    [~, ix] = sortrows([rbin, xang]);

    loopix = find(diff(xang(ix)) < 0);

    cache = struct();
    cache.sz     = sz;
    cache.X      = Xw;
    cache.Y      = Yw;
    cache.camxy  = camw(1:2);
    cache.dx     = dx;
    cache.dy     = dy;
    cache.Xn     = Xn;
    cache.Yn     = Yn;
    cache.rxy    = rxy;
    cache.xang   = xang;
    cache.ix     = ix;
    cache.loopix = loopix;

    if numel(loopix) >= 2
        cache.starts = loopix(1:end-1) + 1;
        cache.stops  = loopix(2:end);
    else
        cache.starts = [];
        cache.stops  = [];
    end
end

rxy    = cache.rxy;
xang   = cache.xang;
ix     = cache.ix;
starts = cache.starts;
stops  = cache.stops;

% ---------------------------------------------------------------------
% Height-dependent part
% ---------------------------------------------------------------------
Zv = Zw(:) - camw(3);
d  = hypot(rxy, Zv);
y  = Zv ./ d;

n   = numel(Zv);
vis = true(n,1);

maxd = max(d,[],'omitnan');
N    = ceil(2*pi / (dx / maxd));

if useSingle
    voxx = single((0:N)' / N);
    voxy = -inf(size(voxx), 'single');
else
    voxx = (0:N)' / N;
    voxy = -inf(size(voxx));
end

% ---------------------------------------------------------------------
% Horizon sweep
% ---------------------------------------------------------------------
for k = 1:numel(starts)
    lp0 = ix(starts(k):stops(k));
    if isempty(lp0)
        continue
    end

    m = numel(lp0);
    lp = zeros(m+2, 1, 'like', lp0);
    lp(2:end-1) = lp0;
    lp(1)       = lp0(end);
    lp(end)     = lp0(1);

    yy = y(lp);
    xx = xang(lp);

    xx(1)   = xx(1) - 1;
    xx(end) = xx(end) + 1;

    % -------------------------------------------------------------
    % Query current horizon at ray angles
    % voxx is already unique/monotonic
    % -------------------------------------------------------------
    Fh = griddedInterpolant(voxx, voxy, 'linear', 'nearest');
    vis(lp0) = Fh(xx(2:end-1)) < yy(2:end-1);

    % -------------------------------------------------------------
    % Update horizon from this ring
    % griddedInterpolant requires unique sample points
    % For duplicate xx, keep maximum yy to preserve upper envelope
    % -------------------------------------------------------------
    [xu, ~, ic] = unique(xx, 'sorted');
    yu = accumarray(ic, yy, [], @max);

    if numel(xu) >= 2
        Fr = griddedInterpolant(xu, yu, 'linear', 'nearest');
        voxy = max(voxy, Fr(voxx));
    end
end

vis = reshape(vis, sz);
end