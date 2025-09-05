function eps_ice = ice_permittivity_tiuri84(frequency, temperature)
%     """Calculates the complex ice dielectric constant depending on the frequency and temperature
%
%     Based on Tiuri et al. (1984). The Complex Dielectric Constant of Snow at Microwave Frequencies.
%     IEEE Journal of Oceanic Engineering, vol. 9, no. 5., pp. 377-382
%
%     :param frequency: frequency in Hz
%     :param temperature: temperature in K
%     :returns: complex permittivity of pure ice
%
%     **Usage example**::
%
%         from smrt.permittivity.ice import ice_permittivity_tiuri84
%         eps_ice = ice_permittivity_tiuri84(frequency=1.9e9, temperature=250)
%
% """

FREEZING_POINT = 273.15;
DENSITY_OF_ICE = 916.7;
tempC = temperature - FREEZING_POINT;

if any(tempC > 0)
    error("The ice temperature must be lower or equal to {FREEZING_POINT}K")
end

% # Units conversion
density_gm3 = DENSITY_OF_ICE * 1e-3;

% # Eq (1) - Real part
Ereal = 1 + 1.7 * density_gm3 + 0.7 * density_gm3^2;

% # Eq (6) - Imaginary part
Eimag = 1.59e6 * ...
    (0.52 * density_gm3 + 0.62*density_gm3^2) * ...
    (frequency^-1 + 1.23e-14 * frequency^.5) * exp(0.036 * tempC);

eps_ice = Ereal + 1j * Eimag;


end