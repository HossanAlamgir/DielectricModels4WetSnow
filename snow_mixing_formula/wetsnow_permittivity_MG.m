function eps_eff = wetsnow_permittivity_MG(freqGHz, density, lwc, ice_permittivity_model, water_permittivity_model)

% % % 3 phase wet snow mixing according to Maxwell Garnett mixing rule as given in Sihvola 1999 book
% % % Alamgir Hossan, JPL, 3/6/2025

%Constatnts
DENSITY_OF_ICE = 916.7;
DENSITY_OF_WATER = 1000;
FREEZING_POINT = 273.15;

if ~exist('temperature','var')
    temperature = FREEZING_POINT;
end

if temperature < FREEZING_POINT && any(lwc > 0)
    error('Liquid water is positive and temperature is set below freezing. This seems incompatible.');
end

% Set default models for water and ice permittivity if not provided
if nargin < 6 || isempty(ice_permittivity_model)
    ice_permittivity_model = @ice_permittivity_maetzler06; % Use a default model
end
if nargin < 5 || isempty(water_permittivity_model)
    water_permittivity_model = @water_permittivity_maetzler87; % Use a default model
end

[frac_volume, fi, fw] = compute_frac_volumes(density, lwc);

eps0=1; % air
eps1=ice_permittivity_model(freqGHz, temperature); %ice
epsd = MG_mixing_first_order(eps0, eps1,fi); % dry inc inc in air env
ew = water_permittivity_model(freqGHz, temperature); % liquid water
eps_eff = MG_mixing_first_order(epsd, ew,fw); % three phase mixing, liq water inc in dry snow env
    

end
