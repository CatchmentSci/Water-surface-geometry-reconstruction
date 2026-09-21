function kh0 = expected_kh(Fr,vel_profile)

g=9.81;


k0init = 1;

switch vel_profile
    case 'constant'
        fun = @(kh)  tanh(kh)/(kh);
end

for i = 1:length(Fr)
        if Fr(i)<=1 % only consider sub-critical Froude number conditions
        kh0(i) = fzero(@(kh) fun(kh) - Fr(i)^2, k0init);
        else
            kh0(i) = NaN;
        end
end


end
