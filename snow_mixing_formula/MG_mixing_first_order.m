function [eps_eff2] = MG_mixing_first_order(epse,epsi,f)

% Maxwell-Garnett Formula for Effective Dielectric Permittivity
% Coder:Alamgir Hossan, JPL, 3/5/2025


% Calculate the effective dielectric permittivity using the Maxwell-Garnett formula
eps_eff1 = epse * ( (epsi + 2*epse - 2*f*(epse - epsi)) ...
    / (epsi + 2*epse + f*(epse - epsi)) ); % Verify with the ref

% Calculate the effective dielectric permittivity using the Maxwell-Garnett
% formula (Sihvola 1999 eq. 3.27)
eps_eff2 = epse + (3*f * epse) * (epsi - epse) / (epsi + 2*epse - f*(epsi - epse));

% the Maxwell-Garnett formula for dioute mixtures (Sihvola 1999 eq. 3.31)
eps_eff3 = epse + (3*f * epse) * (epsi - epse) / (epsi + 2*epse);
end

