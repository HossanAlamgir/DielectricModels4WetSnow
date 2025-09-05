function eps_eff = wetsnow_permittivity_colbeck80_caseII(frequency,density, liquid_water, ice_permittivity_model, water_permittivity_model)


% """Computes the effective permittivity proposed by Colbeck, 1980 for the funicular regime and low dry snow density.
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
if nargin < 5 || isempty(water_permittivity_model)
    water_permittivity_model = @water_permittivity_maetzler87; % Use a default model
end
if nargin < 6 || isempty(ice_permittivity_model)
    ice_permittivity_model = @ice_permittivity_maetzler06; % Use a default model
end


[frac_volume, fi, fw] = compute_frac_volumes(density, liquid_water);

eps1=ice_permittivity_model(frequency, temperature);
eps0=water_permittivity_model(frequency, temperature);
eps2 = 1;

eps_eff=polder_van_santen_three_spherical_components(fi,1-frac_volume,eps0, eps1, eps2);

end
