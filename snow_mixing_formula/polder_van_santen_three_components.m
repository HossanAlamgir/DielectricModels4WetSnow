function eps_eff=polder_van_santen_three_components(f1, f2, eps0, eps1, eps2, A1, A2)


%     """Calculates effective permittivity using Polder and van Santen with three components
%
%     :param f1: fractional volume of component 1
%     :param f2: fractional volume of component 2
%     :param eps0: permittivity of material 0
%     :param eps1: permittivity of material 1
%     :param eps2: permittivity of material 2
%     :param A1: depolarization factor for material 1
%     :param A2: depolarization factor for material 2
%
% """


   if numel(f1) > 1 || numel(f2) > 1
        % If f1 or f2 is an array, vectorize the function
        eps_eff = arrayfun(@(f1_val, f2_val) polder_van_santen_three_components(f1_val, f2_val, eps0, eps1, eps2, A1, A2), f1, f2);
        return;
   end

   % rough first guess
   f0 = 1 - f1 - f2;
   eps_eff0 = f1 * eps1 + f2 * eps2 + f0 * eps0;

    % Define the fn
    pvs = @(x) pvs_equation(x, f1, f2, eps0, eps1, eps2, A1, A2);

    % Solve 
    options = optimset('Display', 'off');  
    [res, ~] = fsolve(pvs, [real(eps_eff0), imag(eps_eff0)], options);

    % Construct the complex permittivity
    eps_eff = complex(res(1), res(2));

end
