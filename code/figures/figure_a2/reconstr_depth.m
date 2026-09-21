function depth = reconstr_depth(wavelength,vel,alpha,d_bounds,vel_profile)

g=9.81;

% initialise depth
depth = nan(length(wavelength),2); % allow for the presence of two possible depth solutions

% wavenumner (rad/m)
wavenumber = 2*pi./wavelength;

if length(vel) < length(wavelength)
    vel = mean(vel)*ones(size(wavelength));
end

m = 2*(1-alpha);
n = (1 - alpha)/alpha;
s = sign(0.5 - n);
nu_num = s*(0.5 - n); % order of the modified Bessel function
nu_den = - s*(0.5 + n); % order of the modified Bessel function


for ik = 1:length(wavenumber)
    k = wavenumber(ik);
    v = vel(ik);

switch vel_profile

    case 'constant'
        % Us = sqrt( (1 + B)/B *g*d *tanh(k*d)/(k*d) ) e.g., Dolcetti et al. (2016), eq. (4)
        fun = @(d)  (g*d)  *tanh(k*d)/(k*d) - v^2;

    case 'linear'
        % (1 - 2\beta)*Us^2 = c_i^2 e.g., Dolcetti et al. (2016), eq. (6)
        fun = @(d) (g*d) /(1-2*((m/2) *tanh(k.*d) ./ (k.*d)))   *tanh(k*d)/(k*d) - v^2;
        %             gd           (1 + B)                     B                             1/(1-2\beta)

    case 'power'
        % k Us^2 = g ( I_{(0.5-n)s}(k*d) / I_{-(0.5-n)s}(k*d) )  e.g., Dolcetti and Garcìa Nava (2019), Closure, eq. (C5)
        fun = @(d) (g/k) * besseli(nu_num, k*d, 1)/besseli(nu_den, k*d, 1) - v^2;

end

d_rec = roots(chebfun(fun, d_bounds));
depth(ik,1:length(d_rec)) = d_rec;

end

depth = depth(:,1)';


end
