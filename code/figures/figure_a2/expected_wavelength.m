function wavelength = expected_wavelength(depth,vel,alpha,vel_profile)

g=9.81;
surf_tens=72.75e-03; %(N/m), ref. 20deg, from Vargaftik et al., 1983 - International tables of the surface tension of water, J.Phys.Chem.Ref.Data
dens=998.2; %(kg/m^3), ref. 20deg

wavelength=nan(size(depth));

k0init = 2*pi;
rootOptions = optimset('Display','off');

switch vel_profile
    case 'constant'
        % Us = sqrt( (1 + B)/B *g*d *tanh(k*d)/(k*d) ) e.g., Dolcetti et al. (2016), eq. (4)
        fun = @(k,d)  (g*d) *(1 + dens*g/(surf_tens*k^2)) / (dens*g/(surf_tens*k^2))  *tanh(k*d)/(k*d);
        %              gd           (1 + B)                     B

    case 'linear'
        m = 2*(1-alpha);
        % (1 - 2\beta)*Us^2 = c_i^2 e.g., Dolcetti et al. (2016), eq. (6)
        fun = @(k,d) (g*d) *(1 + dens*g/(surf_tens*k^2)) / (dens*g/(surf_tens*k^2))   /(1-2*((m/2) *tanh(k.*d) ./ (k.*d)))   *tanh(k*d)/(k*d);
        %             gd           (1 + B)                     B                             1/(1-2\beta)
    case 'power'
        % k Us^2 = g ( I_{(0.5-n)s}(k*d) / I_{-(0.5-n)s}(k*d) )  e.g., Dolcetti and Garcìa Nava (2019), Closure, eq. (C5)
        n = (1 - alpha)/alpha;
        s = sign(0.5 - n);
        nu_num = s*(0.5 - n); % order of the modified Bessel function
        nu_den = - s*(0.5 + n); % order of the modified Bessel function

        fun = @(k,d) (g/k) * besseli(nu_num, k*d, 1)/besseli(nu_den, k*d, 1);
end

for idepth = 1:size(depth,2)
    for ivel = 1:size(vel,1)
        d = depth(ivel,idepth);
        v = vel(ivel,idepth);
        if v^2./(g*d)<=1 % only consider sub-critical Froude number conditions
        k0(ivel,idepth) = fzero(@(k) fun(k,d) - v^2, k0init, rootOptions);
        else
            k0(ivel,idepth) = NaN;
        end
    end
end

wavelength=2*pi./k0;

end
