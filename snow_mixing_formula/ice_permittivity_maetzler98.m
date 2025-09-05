function eps_ice = ice_permittivity_maetzler98(freqGHz, temperature)
%
% """Computes permittivity of ice (accounting for ionic impurities in ice?), equations from Hufford (1991) as given in Maetzler (1998): 'Microwave properties of ice and snow', in B. Schmitt et al. (eds.): 'Solar system ices', p. 241-257, Kluwer.
%
% :param temperature: ice temperature in K
% :param frequency: Frequency in Hz"""

FREEZING_POINT = 273.15;
DENSITY_OF_ICE = 916.7;
tempC = temperature - FREEZING_POINT;

if any(tempC > 0)
    error("The ice temperature must be lower or equal to {FREEZING_POINT}K")
end

epi = 3.1884 + 9.1e-4 * tempC;

% # The Hufford model for the imaginary part:
theta = 300. / temperature - 1;
alpha = (0.00504 + 0.0062 * theta) * np.exp(-22.1 * theta);
beta = (0.502 - 0.131 * theta / (1 + theta)) * 1e-4 + ...
    (0.542e-6 * ((1 + theta) / (theta + 0.0073))^2);

epii = (alpha / freqGHz) + (beta * freqGHz);
eps_ice = epi + epii * 1j;




end