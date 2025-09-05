function eps_eff = wetsnow_permittivity_wiesmann99(freqGHz, density, liquid_water, ice_permittivity_model)


% """Computes the effective permittivity proposed by Colbeck, 1980 for the low porosity.
%
% Colbeck, S. C. (1980). Liquid distribution and the dielectric constant of wet snow.
% Goddard Space Flight Center Microwave Remote Sensing of Snowpack Properties, 21–40.
%
% """


%Constatnts
DENSITY_OF_ICE = 916.7;
DENSITY_OF_WATER = 1000;
FREEZING_POINT = 273.15;

if ~exist('temperature','var')
    temperature = FREEZING_POINT;
end

if temperature < FREEZING_POINT && any(liquid_water > 0)
    error('Liquid water is positive and temperature is set below freezing. This seems incompatible.');
end

% Set default models for water and ice permittivity if not provided
if nargin < 6 || isempty(ice_permittivity_model)
    ice_permittivity_model = @ice_permittivity_maetzler06; % Use a default model
end

[frac_volume, fi, Wi] = compute_frac_volumes(density, liquid_water);

eps0=1;
eps1=ice_permittivity_model(freqGHz, temperature);
eps_dry = polder_van_santen(fi, eps0, eps1);

Aa = 0.005;    % depolarisation factors of prolate
Ab = 0.4975;   % water inclusion (Matzler 1987)
Ac = Ab;

eps_sw = 88;
eps_inf_w = 4.9;
f0w = 9e9;  % GHz

eps_eff = 0;
frequency = freqGHz * 1e9;
for Ak = [Aa, Ab, Ac]
    eps_s_k = Wi / 3 * (eps_sw - eps_dry) / (1 + Ak * (eps_sw / eps_dry - 1));
    eps_inf_k = Wi / 3 * (eps_inf_w - eps_dry) / (1 + Ak * (eps_inf_w / eps_dry - 1));
    f0_k = f0w * (1 + Ak * (eps_sw - eps_inf_w) / (eps_dry + Ak * (eps_inf_w - eps_dry)));

    eps_k = eps_inf_k + (eps_s_k - eps_inf_k) / (1 - 1j * frequency / f0_k);

    eps_eff = eps_eff + eps_k;

end
eps_eff = eps_dry + eps_eff;

end
