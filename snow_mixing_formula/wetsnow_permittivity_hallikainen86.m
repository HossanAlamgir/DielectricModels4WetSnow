function eps = wetsnow_permittivity_hallikainen86( freqGHz,dry_snow_density_gcm3, mv)
%
% """Computes the effective permittivity of a snow mixture calculated with the Modified Debye model by Hallikainen 1986
%
% The implemented equation are 10, 11 and 13a-c.
%
% The validity of the model is: frequency between 3 and 37GHz;
% mv between 1% and 12%;
% dry_snow_density between 0.09 and 0.38g/cm3.
%
% The implementation of this function follows the equations formulation of the original paper
% Hallikainen, M., F. Ulaby, and M. Abdelrazik, “Dielectric properties of snow in 3 to 37 GHz range,”
% IEEE Trans. on Antennasand Propagation,Vol. 34, No. 11, 1329–1340, 1986. DOI: 10.1109/TAP.1986.1143757
% Anyway this formulation does not allow the reproduction of the results as reported in the paper.
% A new formulation of eq. 12a have been presented in the book
% Microwave Radar and Radiometric Remote Sensing by Ulaby et al. 2014 from which the SMRT function
% wetsnow_permittivity_hallikainen86_ulaby14 have been implemented. The users are pointed to that definition.
%
% smrt_warn("This model cannot reproduce the results of the original paper. You may want to use wetsnow_permittivity_hallikainen86_ulaby14 instead.")
% """

%Constatnts
DENSITY_OF_ICE = 916.7;
DENSITY_OF_WATER = 1000;
FREEZING_POINT = 273.15;

% [frac_volume, fi, fw] = compute_frac_volumes(density, liquid_water);
% 
% % # fractional volume of water in %
% mv = 100 * fw;
% 
% % # Eq 3 in H86 defines the dry snow by (snow density - mass of water per volume of snow) / (1 - volume fo water per volume of snow)
% dry_snow_density_gcm3 = 1e-3 * (density - DENSITY_OF_WATER * fw) / (1 - fw);


% # Eq 13
A1 = 0.78 + 0.03 * freqGHz - 0.58e-3 * freqGHz^2;
A2 = 0.97 - 0.39e-2 * freqGHz + 0.39e-3 * freqGHz^2;
B1 = 0.31 - 0.05 * freqGHz + 0.87e-3 * freqGHz^2;

% # Eq 12 (different from Ulaby14)
A = 1 + 1.83 * dry_snow_density_gcm3 + 0.02 * A1 * mv^1.015 + B1;
B = 0.073 * A1;
C = 0.073 * A2;
x = 1.31;

freq0 = 9.07;%  # GHz

eps_ws_r = A + B * mv^x / (1 + (freqGHz / freq0)^2);

eps_ws_i = C * mv^x * (freqGHz / freq0) / (1 + (freqGHz / freq0)^2);

eps = eps_ws_r + 1j * eps_ws_i;


end
