classdef camera

    properties
        xyz = [0 0 0];                  % world coordinates of camera
        imgsz = [100 100];              % [rows, cols]
        viewdir = [0 0 0];              % [yaw pitch roll]
        f = [5000 5000];                % [fx fy]
        c = [50 50];                    % [cx cy]
        k = [0 0 0 0 0 0];              % [k1..k6]
        p = [0 0];                      % [p1 p2]
    end

    properties (Dependent)
        R
        fullmodel
    end

    methods
        function cam = camera(varargin)
            % Constructor
            if nargin == 0
                return
            end

            if nargin == 1
                cam.fullmodel = varargin{1};
                return
            end

            if nargin < 7
                varargin{7} = [];
            end

            cam.xyz     = varargin{1};
            cam.imgsz   = varargin{2};
            cam.viewdir = varargin{3};
            cam.f       = varargin{4};
            cam.c       = varargin{5};
            cam.k       = varargin{6};
            cam.p       = varargin{7};

            if numel(cam.imgsz) < 2
                error('malformed image size.');
            end

            cam.imgsz(3:end) = [];
            cam.f(end+1:2) = cam.f(end);

            if isempty(cam.c)
                cam.c = (cam.imgsz([2 1]) + 1) / 2;
            end

            cam.k(end+1:6) = 0;
            cam.p(end+1:2) = 0;
        end

        function value = get.R(cam)
            value = camera.rotationFromViewdir(cam.viewdir);
        end

        function value = get.fullmodel(cam)
            value = [cam.xyz, cam.imgsz, cam.viewdir, cam.f, cam.c, cam.k, cam.p];
        end

        function cam = set.fullmodel(cam, value)
            cam.xyz     = value(1:3);
            cam.imgsz   = value(4:5);
            cam.viewdir = value(6:8);
            cam.f       = value(9:10);
            cam.c       = value(11:12);
            cam.k       = value(13:18);
            cam.p       = value(19:20);
        end

        function [uv, depth, inframe] = project(cam, xyz)
            % Project world coordinates to image coordinates

            if size(xyz,2) ~= 3 && size(xyz,1) == 3
                xyz = xyz.';
            end

            [uv, depth, inframe] = camera.projectNumeric( ...
                xyz, ...
                cam.xyz, ...
                cam.R, ...
                cam.f, ...
                cam.c, ...
                cam.k, ...
                cam.p, ...
                cam.imgsz ...
            );
        end

        function prep = prepareDEMInverse(cam, X, Y, Z)
            % Prepare static DEM-grid state for repeated inverse projections.

            if nargin < 4
                Z = [];
            end

            if isempty(X) || isempty(Y)
                error('X and Y must be supplied.');
            end

            [Xo, Yo, flipud_flag, fliplr_flag] = camera.orientGrid(X, Y);

            prep = struct;
            prep.prepType = 'cameraDEMPrep';
            prep.X = Xo;
            prep.Y = Yo;
            prep.flipud_flag = flipud_flag;
            prep.fliplr_flag = fliplr_flag;
            prep.xvec = Xo(1,:);
            prep.yvec = Yo(:,1);
            prep.hasGridInterpolant = exist('griddedInterpolant', 'class') ~= 0;
            prep.hasScatteredInterpolant = exist('scatteredInterpolant', 'file') > 1;
            prep.Z = [];
            prep.visible = [];
            prep.zfun = [];

            if nargin >= 4 && ~isempty(Z)
                prep = cam.updatePreparedDEM(prep, Z);
            end
        end

        function prep = updatePreparedDEMInterpOnly(cam, prep, Z) %#ok<INUSD>
            % Fast DEM update for trial evaluations:
            % reuses existing visibility and updates only DEM/interpolant values.
            %
            % This is intended for provisional row tests where the DEM perturbation
            % is tiny and a full voxelviewshed recomputation is too expensive.
            %
            % IMPORTANT:
            %   prep.visible must already come from a full updatePreparedDEM call.

            if ~isstruct(prep) || ~isfield(prep, 'prepType') || ~strcmp(prep.prepType, 'cameraDEMPrep')
                error('prep must be a struct produced by prepareDEMInverse.');
            end

            if isempty(prep.visible)
                error('prep.visible is empty. Run updatePreparedDEM first.');
            end

            Z_use = Z;
            if prep.flipud_flag
                Z_use = flipud(Z_use);
            end
            if prep.fliplr_flag
                Z_use = fliplr(Z_use);
            end

            % Reuse existing visibility mask from the current accepted DEM state
            Z_vis = Z_use ./ prep.visible;

            prep.Z = Z_vis;

            if prep.hasGridInterpolant
                if isfield(prep, 'zfun') && isa(prep.zfun, 'griddedInterpolant')
                    prep.zfun.Values = Z_vis.';
                else
                    prep.zfun = griddedInterpolant({prep.xvec, prep.yvec}, Z_vis.', 'linear', 'none');
                end
            else
                prep.zfun = [];
            end
        end






        function prep = updatePreparedDEM(cam, prep, Z)
            % Update the DEM values and interpolant/visibility on a prepared grid.

            if ~isstruct(prep) || ~isfield(prep, 'prepType') || ~strcmp(prep.prepType, 'cameraDEMPrep')
                error('prep must be a struct produced by prepareDEMInverse.');
            end

            Z_use = Z;
            if prep.flipud_flag
                Z_use = flipud(Z_use);
            end
            if prep.fliplr_flag
                Z_use = fliplr(Z_use);
            end

            visible = voxelviewshed(prep.X, prep.Y, Z_use, cam.xyz);
            Z_vis = Z_use ./ visible;

            prep.Z = Z_vis;
            prep.visible = visible;

            if prep.hasGridInterpolant
                if isfield(prep, 'zfun') && isa(prep.zfun, 'griddedInterpolant')
                    prep.zfun.Values = Z_vis.';
                else
                    prep.zfun = griddedInterpolant({prep.xvec, prep.yvec}, Z_vis.', 'linear', 'none');
                end
            else
                prep.zfun = [];
            end
        end

        function prep = perturbPreparedDEM(cam, prep, iy, ix, dz)
            % Copy a prepared DEM and apply a single-cell perturbation safely.

            if isempty(prep.Z)
                error('prep does not contain DEM values. Call updatePreparedDEM first.');
            end

            Z_new = prep.Z;
            Z_new(iy, ix) = Z_new(iy, ix) + dz;
            prep = cam.updatePreparedDEM(prep, Z_new);
        end

        function xyz = invproject(cam, uv, X, Y, Z, xy0)
            % Inverse projection from image to world coordinates.
            %
            % Supported call patterns:
            %   xyz = cam.invproject(uv)
            %   xyz = cam.invproject(uv, X, Y, Z)
            %   xyz = cam.invproject(uv, X, Y, Z, xy0)
            %   xyz = cam.invproject(uv, prep)
            %   xyz = cam.invproject(uv, prep, xy0)

            nanix = any(isnan(uv), 2);
            anynans = any(nanix);

            if anynans
                uv_valid = uv(~nanix, :);
            else
                uv_valid = uv;
            end

            nUV = size(uv_valid, 1);

            cam_xyz = cam.xyz;
            R = cam.R;
            f = cam.f;
            c = cam.c;
            k = cam.k;
            p = cam.p;
            imgsz = cam.imgsz;
            hasDist = any(k ~= 0) || any(p ~= 0);

            if nargin == 2
                % Exact non-distorted inverse, then refine for distortion if needed

                depth0 = 1000;
                xyz_valid = [(uv_valid(:,1)-c(1))/f(1), ...
                             (uv_valid(:,2)-c(2))/f(2)] * depth0;
                xyz_valid(:,3) = depth0;
                xyz_valid = xyz_valid * R;
                xyz_valid = xyz_valid + cam_xyz;

                if hasDist
                    E = [1 0 0; 0 1 0] * R;  % perturbation directions

                    for ii = 1:nUV
                        base_xyz = xyz_valid(ii,:);
                        uv_i = uv_valid(ii,:);

                        m = LMFnlsq( ...
                            @(m) camera.invprojectResidualRay( ...
                                m, base_xyz, E, uv_i, cam_xyz, R, f, c, k, p, imgsz), ...
                            [0; 0]);

                        xyz_valid(ii,:) = base_xyz + m(:).' * E;
                    end
                end

                xyz = camera.reinsertNaNs(xyz_valid, nanix, size(uv,1), anynans);
                return
            end

            if nargin >= 3 && isstruct(X) && isfield(X, 'prepType') && strcmp(X.prepType, 'cameraDEMPrep')
                prep = X;
                if nargin >= 4
                    xy0 = Y;
                else
                    xy0 = [];
                end
                xyz_valid = camera.invprojectPrepared( ...
                    uv_valid, prep, xy0, cam_xyz, R, f, c, k, p, imgsz);
                xyz = camera.reinsertNaNs(xyz_valid, nanix, size(uv,1), anynans);
                return
            end

            if nargin < 5
                error('DEM-constrained inverse projection requires X, Y and Z.');
            end

            prep = cam.prepareDEMInverse(X, Y, Z);

            if nargin < 6
                xy0 = [];
            end

            xyz_valid = camera.invprojectPrepared( ...
                uv_valid, prep, xy0, cam_xyz, R, f, c, k, p, imgsz);

            xyz = camera.reinsertNaNs(xyz_valid, nanix, size(uv,1), anynans);
        end

        function [result, rmse, AIC] = optimizecam(cam, xyz, uv, freeparams)
            % Tune the camera so projected xyz matches uv

            nanrows = any(isnan(xyz),2) | any(isnan(uv),2);
            xyz(nanrows,:) = [];
            uv(nanrows,:) = [];

            fullmodel0 = cam.fullmodel;

            freeparams = ~(freeparams(:)==0 | freeparams(:)=='0')';
            paramix = find(freeparams);
            Nfree = numel(paramix);
            mbest = zeros(Nfree,1);

            if size(uv,2) == 3
                weights = uv(:,3);
                uv_target = uv(:,1:2);
                misfit = @(m) camera.optimizeResidualWeighted(fullmodel0, paramix, m, xyz, uv_target, weights);
            else
                uv_target = uv;
                misfit = @(m) camera.optimizeResidual(fullmodel0, paramix, m, xyz, uv_target);
            end

            r0 = misfit(mbest);
            if any(isnan(r0))
                error('All GCPs must be infront of the initial camera location for optimizecam to work.');
            end

            [mbest, RSS] = LMFnlsq(misfit, mbest);

            Nuv = size(uv_target,1);
            rmse = sqrt(RSS / Nuv);
            AIC = numel(uv_target) * log(RSS / numel(uv_target)) + 2 * Nfree;

            fullmodel_best = fullmodel0;
            fullmodel_best(paramix) = fullmodel_best(paramix) + mbest(:).';
            result = camera(fullmodel_best);
        end
    end

    methods (Static, Access = private)
        function R = rotationFromViewdir(viewdir)
            C = cos(viewdir);
            S = sin(viewdir);

            R = [ ...
                S(3).*S(2).*C(1)-C(3).*S(1),  S(3).*S(2).*S(1)+C(3).*C(1),  S(3).*C(2); ...
                C(3).*S(2).*C(1)+S(3).*S(1),  C(3).*S(2).*S(1)-S(3).*C(1),  C(3).*C(2); ...
                C(2).*C(1),                   C(2).*S(1),                  -S(2) ...
            ];

            R(1:2,:) = -R(1:2,:);
        end

        function restoreLMWarningState(warnStateSingular, warnStateNearly)
            warning(warnStateSingular.state, 'MATLAB:singularMatrix');
            warning(warnStateNearly.state,   'MATLAB:nearlySingularMatrix');
        end

        function [Xo, Yo, flipud_flag, fliplr_flag] = orientGrid(X, Y)
            Xo = X;
            Yo = Y;
            flipud_flag = false;
            fliplr_flag = false;

            if Yo(2,2) < Yo(1,1)
                Xo = flipud(Xo);
                Yo = flipud(Yo);
                flipud_flag = true;
            end

            if Xo(2,2) < Xo(1,1)
                Xo = fliplr(Xo);
                Yo = fliplr(Yo);
                fliplr_flag = true;
            end
        end

        function xyz = reinsertNaNs(xyz_valid, nanix, nRows, anynans)
            if anynans
                xyz = nan(nRows, 3);
                xyz(~nanix,:) = xyz_valid;
            else
                xyz = xyz_valid;
            end
        end

        function xyz_valid = invprojectPrepared(uv_valid, prep, xy0, cam_xyz, R, f, c, k, p, imgsz)
            nUV = size(uv_valid, 1);
            xyz_valid = nan(nUV, 3);

            if isempty(prep.Z)
                error('Prepared DEM state does not contain Z values. Call updatePreparedDEM first.');
            end

            if isempty(xy0)
                pts = [prep.X(prep.visible(:)), prep.Y(prep.visible(:)), prep.Z(prep.visible(:))];
                [uv0, ~, inframe] = camera.projectNumeric( ...
                    pts, cam_xyz, R, f, c, k, p, imgsz);

                uv0 = [uv0, pts];
                uv0 = uv0(inframe, :);

                if prep.hasScatteredInterpolant
                    Xscat = scatteredInterpolant(uv0(:,1), uv0(:,2), uv0(:,3));
                    Yscat = scatteredInterpolant(uv0(:,1), uv0(:,2), uv0(:,4));
                    Zscat = scatteredInterpolant(uv0(:,1), uv0(:,2), uv0(:,5));
                else
                    % fallback for older MATLAB
                    Xscat = TriScatteredInterp(uv0(:,1), uv0(:,2), uv0(:,3)); %#ok<REMFF1>
                    Yscat = TriScatteredInterp(uv0(:,1), uv0(:,2), uv0(:,4)); %#ok<REMFF1>
                    Zscat = TriScatteredInterp(uv0(:,1), uv0(:,2), uv0(:,5)); %#ok<REMFF1>
                end

                xyz_valid(:,1) = Xscat(uv_valid(:,1), uv_valid(:,2));
                xyz_valid(:,2) = Yscat(uv_valid(:,1), uv_valid(:,2));
                xyz_valid(:,3) = Zscat(uv_valid(:,1), uv_valid(:,2));
                return
            end

            if prep.hasGridInterpolant
                zfun = prep.zfun;
            else
                zfun = @(x,y) interp2(prep.X, prep.Y, prep.Z, x, y);
            end

            % Convert singular LM warnings into catchable errors for this solve block.
            warnStateSingular = warning('query', 'MATLAB:singularMatrix');
            warnStateNearly   = warning('query', 'MATLAB:nearlySingularMatrix');
            warning('error', 'MATLAB:singularMatrix');
            warning('error', 'MATLAB:nearlySingularMatrix');
            cleanupWarn = onCleanup(@() camera.restoreLMWarningState( ...
                warnStateSingular, warnStateNearly)); %#ok<NASGU>

            for ii = 1:nUV
                uv_i = uv_valid(ii,1:2);
                xy_init = xy0(ii,1:2);

                % Reject invalid warm starts immediately.
                if any(~isfinite(xy_init))
                    continue
                end

                z_init = zfun(xy_init(1), xy_init(2));
                if ~isfinite(z_init)
                    continue
                end

                r_init = camera.invprojectResidualDEM( ...
                    xy_init(:), zfun, uv_i, cam_xyz, R, f, c, k, p, imgsz);

                if any(~isfinite(r_init))
                    continue
                end

                try
                    xy_sol = LMFnlsq( ...
                        @(xy) camera.invprojectResidualDEM( ...
                        xy, zfun, uv_i, cam_xyz, R, f, c, k, p, imgsz), ...
                        xy_init(:));

                    if numel(xy_sol) < 2 || any(~isfinite(xy_sol(1:2)))
                        continue
                    end

                    z_sol = zfun(xy_sol(1), xy_sol(2));
                    if ~isfinite(z_sol)
                        continue
                    end

                    err = camera.invprojectResidualDEM( ...
                        xy_sol, zfun, uv_i, cam_xyz, R, f, c, k, p, imgsz);

                    if any(~isfinite(err)) || sum(err.^2) > 2^2
                        continue
                    end

                    xyz_valid(ii,1:2) = xy_sol(1:2).';
                    xyz_valid(ii,3)   = z_sol;

                catch
                    % Singular / nearly singular / failed LM solve:
                    % leave this point as NaN and continue.
                end
            end
        end

        function [uv, depth, inframe] = projectNumeric(xyz, cam_xyz, R, f, c, k, p, imgsz)
            % Fast numeric projector used internally

            if size(xyz,2) ~= 3 && size(xyz,1) == 3
                xyz = xyz.';
            end

            xyz = xyz - cam_xyz;
            xyz = xyz * R.';

            depth = xyz(:,3);
            xy = xyz(:,1:2) ./ depth;

            hasDist = any(k ~= 0) || any(p ~= 0);

            if hasDist
                x = xy(:,1);
                y = xy(:,2);

                r2 = x.^2 + y.^2;
                r2(r2 > 4) = 4;

                r4 = r2.^2;
                r6 = r4 .* r2;

                num = 1 + k(1)*r2 + k(2)*r4 + k(3)*r6;

                if any(k(4:6) ~= 0)
                    den = 1 + k(4)*r2 + k(5)*r4 + k(6)*r6;
                    a = num ./ den;
                else
                    a = num;
                end

                xyprod = x .* y;
                x2 = x.^2;
                y2 = y.^2;

                % Brown-Conrady tangential distortion:
                % xd = x + 2*p1*x*y + p2*(r^2 + 2*x^2)
                % yd = y + p1*(r^2 + 2*y^2) + 2*p2*x*y
                xd = a .* x + 2*p(1)*xyprod + p(2) * (r2 + 2*x2);
                yd = a .* y + p(1) * (r2 + 2*y2) + 2*p(2)*xyprod;

                uv = [f(1)*xd + c(1), f(2)*yd + c(2)];
            else
                uv = [f(1)*xy(:,1) + c(1), f(2)*xy(:,2) + c(2)];
            end

            uv(depth <= 0, :) = NaN;

            if nargout > 2
                inframe = (depth > 0) & ...
                          (uv(:,1) >= 1) & (uv(:,2) >= 1) & ...
                          (uv(:,1) <= imgsz(2)) & (uv(:,2) <= imgsz(1));
            end
        end

        function r = invprojectResidualRay(m, base_xyz, E, uv_i, cam_xyz, R, f, c, k, p, imgsz)
            xyzp = base_xyz + m(:).' * E;
            uvp = camera.projectNumeric(xyzp, cam_xyz, R, f, c, k, p, imgsz);
            r = (uvp - uv_i).';
        end

        function r = invprojectResidualDEM(xy, zfun, uv_i, cam_xyz, R, f, c, k, p, imgsz)
            z = zfun(xy(1), xy(2));
            uvp = camera.projectNumeric([xy(1), xy(2), z], cam_xyz, R, f, c, k, p, imgsz);
            r = (uvp - uv_i).';
        end

        function r = optimizeResidual(fullmodel0, paramix, m, xyz, uv)
            fullmodel = fullmodel0;
            fullmodel(paramix) = fullmodel(paramix) + m(:).';

            cam_xyz = fullmodel(1:3);
            imgsz   = fullmodel(4:5);
            viewdir = fullmodel(6:8);
            f       = fullmodel(9:10);
            c       = fullmodel(11:12);
            k       = fullmodel(13:18);
            p       = fullmodel(19:20);

            R = camera.rotationFromViewdir(viewdir);
            uvhat = camera.projectNumeric(xyz, cam_xyz, R, f, c, k, p, imgsz);
            r = reshape(uvhat - uv, [], 1);
        end

        function r = optimizeResidualWeighted(fullmodel0, paramix, m, xyz, uv, w)
            fullmodel = fullmodel0;
            fullmodel(paramix) = fullmodel(paramix) + m(:).';

            cam_xyz = fullmodel(1:3);
            imgsz   = fullmodel(4:5);
            viewdir = fullmodel(6:8);
            f       = fullmodel(9:10);
            c       = fullmodel(11:12);
            k       = fullmodel(13:18);
            p       = fullmodel(19:20);

            R = camera.rotationFromViewdir(viewdir);
            uvhat = camera.projectNumeric(xyz, cam_xyz, R, f, c, k, p, imgsz);
            r = reshape((uvhat - uv) .* w(:,[1 1]), [], 1);
        end
    end
end
