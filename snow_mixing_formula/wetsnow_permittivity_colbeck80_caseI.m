function eps_eff = wetsnow_permittivity_colbeck80_caseI(freqGHz, density, liquid_water, ice_permittivity_model, water_permittivity_model)


%     """Computes the effective permittivity proposed by Colbeck, 1980 for the pendular regime.
%
%     Colbeck, S. C. (1980). Liquid distribution and the dielectric constant of wet snow.
%     Goddard Space Flight Center Microwave Remote Sensing of Snowpack Properties, 21–40.
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


Ac = 0.422;% page 24
Asnow = [(1 - Ac) / 2, (1 - Ac) / 2, 0.422];

n = 3.5;
% # for n = 3.5 (page 4), we read in Fig 2 page 31:
m = 0.072;%  # this value is different from Löwe et al. TC (2013) and I don't know why.
Ac = 1 / (1 + 2 / m);

% [Aa,Ab,Ac] = depolarizarion_factor(n);%Jones and Friedman, 2000

Awater = [(1 - Ac) / 2, (1 - Ac) / 2, Ac];

[frac_volume, fi, fw] = compute_frac_volumes(density, liquid_water);

eps0=1;
eps1=ice_permittivity_model(freqGHz, temperature);
eps2=water_permittivity_model(freqGHz, temperature);

eps_eff=polder_van_santen_three_components(fi,fw,eps0, eps1, eps2, Asnow, Awater);

end
