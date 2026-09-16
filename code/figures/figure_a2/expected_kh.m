function kh0 = expected_kh(Fr,alpha,vel_profile)

g=9.81;
surf_tens=72.75e-03; %(N/m), ref. 20deg, from Vargaftik et al., 1983 - International tables of the surface tension of water, J.Phys.Chem.Ref.Data
dens=998.2; %(kg/m^3), ref. 20deg


k0init = 1;

switch vel_profile
    case 'constant'
        % Us = sqrt( (1 + B)/B *g*d *tanh(k*d)/(k*d) ) e.g., Dolcetti et al. (2016), eq. (4)
        fun = @(kh)  tanh(kh)/(kh);
        %              gd           (1 + B)                     B
end

for i = 1:length(Fr)
        if Fr(i)<=1 % only consider sub-critical Froude number conditions
        kh0(i) = fzero(@(kh) fun(kh) - Fr(i)^2, k0init);
        else 
            kh0(i) = NaN;
        end
end


end