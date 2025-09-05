function eps_ice = ice_permittivity_maetzler06(freqGHz, temperature)

% """ Calculates the complex ice dielectric constant depending on the frequency and temperature
%
% Based on Mätzler, C. (2006). Thermal Microwave Radiation: Applications for Remote Sensing p456-461
% This is the default model used in smrt.inputs.make_medium.make_snow_layer().
%
% :param frequency: frequency in Hz
% :param temperature: temperature in K
% :returns: complex permittivity of pure ice
%
% **Usage example**::
%
% from smrt.permittivity.ice import ice_permittivity_maetzler06
% eps_ice = ice_permittivity_maetzler06(frequency=18e9, temperature=270)
%
% .. note::
%
% Ice permittivity is automatically calculated in smrt.inputs.make_medium.make_snow_layer() and
% is not set by the electromagnetic model module. An alternative
% to ice_permittivity_maetzler06 may be specified as an argument to the make_snow_layer
%     function. The usage example is provided for external reference or testing purposes.
%
%         """

FREEZING_POINT = 273.15;
DENSITY_OF_ICE = 916.7;
tempC = temperature - FREEZING_POINT;

if any(tempC > 0)
    error("The ice temperature must be lower or equal to {FREEZING_POINT}K")
end

Ereal = 3.1884 + 9.1e-4 * tempC;

theta = 300.0 / temperature - 1.0;
alpha = (0.00504 + 0.0062 * theta) * exp(-22.1 * theta);

B1 = 0.0207;
B2 = 1.16e-11;
b = 335;
deltabeta = exp(- 9.963 + 0.0372 * tempC);
betam = (B1 / temperature) * (exp(b / temperature) / ((exp(b / temperature) - 1)^2)) + B2 * freqGHz^2;
beta = betam + deltabeta;

Eimag = alpha / freqGHz + beta * freqGHz;

eps_ice = Ereal + 1j * Eimag;


end