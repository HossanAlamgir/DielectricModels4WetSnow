function residual = pvs_equation_spherical(x, f1, f2, eps0, eps1, eps2)
% Polder-Van Santen equation
eps_eff = complex(x(1), x(2));
residual = eps_eff * (1 - 3 * f2 * (eps2 - eps0) / (2 * eps_eff + eps2) - 3 * f1 * (eps1 - eps0) / (2 * eps_eff + eps1)) - eps0;
% Return real and imaginary parts of the residual
residual = [real(residual), imag(residual)];

end
