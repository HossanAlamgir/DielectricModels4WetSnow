function eps_eff = wetsnow_permittivity_Matzler(freqGHz, density, liquid_water, ice_permittivity_model, water_permittivity_model)

% """Computes the effective permittivity of a snow mixture as calculated in MEMLS using Maxwell-Garnett Mixing rule of water in dry snow
% for prolate spheroidal water with experimentally determined. Dry snow permittivity is here determined with Polder van Santen.
% 
%     """
%     # %   depolarisation factors.
%     # %   calculates complex dielectric constant of wet snow
%     # %   using Maxwell-Garnett Mixing rule of water in dry snow
%     # %   for prolate spheroidal water with experimentally determined
%     # %   depolarisation factors.
%     # %   Water temperature is at 273.15 K, with epsilon
%     # %   of water from Liebe et al. 1991.
%     # %       epsd:  complex epsilon of dry snow
%     # %       f:   frequency [GHz]
%     # %       Ti:  physical snow temperature [K]
%     # %       Wi:  wetness [volume fraction]
%     # %
%     # %   Version history:
%     # %      1.0    wi 15.7.95
%     # %      2.0    ma 31.5.2005: Wi is volume fraction (not %)
%     # %      3.0    ma 2.4.2007 : adjustments, new function name
%     # %   Uses: epswater (since Version 3)
%     # %
%     # %   Copyright (c) 1997 by the Institute of Applied Physics,
%     # %   University of Bern, Switzerland
%     # %


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

[frac_volume, fi, Wi] = compute_frac_volumes(density, liquid_water);

eps0=1;
eps1=ice_permittivity_model(freqGHz, temperature);
% epsd = polder_van_santen(fi, eps0, eps1);
TdryC = -0.5;
[epsd] = diel_dry_snow_matzler_06(TdryC, density/1000, freqGHz );
ew = water_permittivity_model(freqGHz, temperature);

Aa = 0.005;    % depolarisation factors of prolate
Ab = 0.4975;   % water inclusion (Matzler 1987)
% Ac = Ab;


Ka = epsd / (epsd + Aa * (ew - epsd));
Kb = epsd / (epsd + Ab * (ew - epsd));
K = (Ka + 2 * Kb) / 3;
epsz = (1 - Wi) * epsd + Wi * ew * K;
epsn = 1 - Wi * (1 - K);
eps_eff = epsz / epsn;  % Maxwell-Garnett Mixing of water in dry snow
    

end
