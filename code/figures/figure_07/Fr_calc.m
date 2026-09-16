function Fr = Fr_calc(kd, velProfile, alpha)
%FR_CALC Calculate Froude number for a specified velocity profile.

switch velProfile

    case 'constant'

        Fr = sqrt(tanh(kd) ./ kd);

    case 'linear'

        m = 2 * (1 - alpha);
        Fr = sqrt(tanh(kd) ./ ...
            (kd - m * tanh(kd)));

    case 'power'

        n = (1 - alpha) / alpha;
        s = sign(0.5 - n);

        nuNum = s * (0.5 - n);
        nuDen = -s * (0.5 + n);

        Fr = sqrt( ...
            besseli(nuNum, kd, 1) ./ ...
            besseli(nuDen, kd, 1) ./ kd);

    otherwise
        error('Unknown velocity profile: %s', velProfile);

end

end