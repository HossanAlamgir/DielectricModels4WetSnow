function eps_eff = wetsnow_permittivity_three_component_polder_van_santen(freqGHz, density, liquid_water, ice_permittivity_model, water_permittivity_model)

% """Computes the effective permittivity of a snow mixture using the three components polder_van_santen, assuming spherical inclusions
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
if nargin < 5 || isempty(water_permittivity_model)
    water_permittivity_model = @water_permittivity_maetzler87; % Use a default model
end

[frac_volume, fi, fw] = compute_frac_volumes(density, liquid_water);

eps0 = 1;
eps1=ice_permittivity_model(freqGHz, temperature);
eps2=water_permittivity_model(freqGHz, temperature);


eps_eff=polder_van_santen_three_spherical_components(fi,fw,eps0, eps1, eps2);


end
