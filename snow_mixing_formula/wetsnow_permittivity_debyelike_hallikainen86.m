function eps = wetsnow_permittivity_debyelike_hallikainen86(freqGHz, dry_snow_density_gcm3, mv)
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
% below 15 GHz, we may set A1 = A2 = 1 .O and Bl = 0 for Debye-like model


% # Eq 13 Discussion
A1 = 1;
A2 = 1;
B1 = 0;

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
