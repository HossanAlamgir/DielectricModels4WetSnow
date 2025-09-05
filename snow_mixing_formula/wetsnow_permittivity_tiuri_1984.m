function eps_eff = wetsnow_permittivity_tiuri_1984(f,Ps, mv)

% Tiuri et al. (1984) Empirical Model
% ϵ_s^'=1+1.7 ρ_d+0.7 〖ρ_d〗^2+8.7 W+70 W^2                                  
% ϵ_s''=f/10^9 (0.9 W+7.5 W^2),  f=500-1000 MHz
%Input Variables:
%f: frequency in GHz
%Ps: Dry Snow Density (g/cm^3)
%mv: volumetric water content (0<mv<6%)
%Output Products:
%eps_eff:complex dielectric constant
% Alamgir Hossan, JPL, 3/5/2025


epsr = 1+1.7*Ps+0.7*Ps^2+8.7*mv+70*mv^2;
epsi = f*(0.9*mv+7.5*mv^2) ;

eps_eff = complex(epsr,epsi);
end