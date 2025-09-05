function residual = pvs_equation(x, f1, f2, eps0, eps1, eps2, A1, A2)
% Polder-Van Santen equation
eps_eff = complex(x(1), x(2));
sum_A2 = sum(1 ./ (eps_eff + A2 .* (eps2 - eps_eff)));
sum_A1 = sum(1 ./ (eps_eff + A1 .* (eps1 - eps_eff)));

residual_complex = eps_eff * (1 - (1 / 3) * f2 * (eps2 - eps0) * sum_A2 - (1 / 3) * f1 * (eps1 - eps0) * sum_A1) - eps0;

% Return real and imaginary parts of the residual
residual = [real(residual_complex), imag(residual_complex)];

end
