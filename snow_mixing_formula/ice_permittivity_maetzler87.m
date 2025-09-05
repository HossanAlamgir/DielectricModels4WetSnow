function eps_ice = ice_permittivity_maetzler87(freqGHz, temperature)

% """Calculates the complex ice dielectric constant depending on the frequency and temperature
% 
% Based on Mätzler, C. and Wegmüller (1987). Dielectric properties of fresh-water ice at microwave frequencies.
% J. Phys. D: Appl. Phys. 20 (1987) 1623-1630.
% 
% :param frequency: frequency in Hz
% :param temperature: temperature in K
% :returns: complex permittivity of pure ice
% 
% **Usage example**::
% 
% from smrt.permittivity.ice import ice_permittivity_maetzler87
% eps_ice = ice_permittivity_maetzler87(frequency=18e9, temperature=270)
% 
% .. note::
% 
% This is only suitable for testing at -5 deg C and -15 deg C. If used at other temperatures
% a warning will be displayed.
% 
% """

FREEZING_POINT = 273.15;
DENSITY_OF_ICE = 916.7;
tempC = temperature - FREEZING_POINT;

if any(tempC > 0)
    error("The ice temperature must be lower or equal to {FREEZING_POINT}K")
end

    % # Equation 10
    Ereal = 3.1884 + 9.1e-4 * tempC;

    if (temperature - FREEZING_POINT) < -10
        A = 3.5e-4;
        B = 3.6e-5;
        C = 1.2;
    else
        A = 6e-4;
        B = 6.5e-5;
        C = 1.07;
    end
    % # Equation 13
    Eimag = A / freqGHz + B * freqGHz^C;
    % # Issue warning if temperature different from values in paper
    if ~ismember(temperature, [FREEZING_POINT - 5, FREEZING_POINT - 15])
        warning('Strictly, this permittivity formulation was proposed for -5 and -15 deg C. It is recommended to use another formulation if this is not for testing purpose.');
    end


eps_ice = Ereal + 1j * Eimag;


end