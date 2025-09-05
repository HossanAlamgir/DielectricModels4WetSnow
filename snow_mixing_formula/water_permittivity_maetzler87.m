function Ew = water_permittivity_maetzler87(freqGHz, temperature)

% """Calculates the complex water dielectric constant depending on the frequency and temperature
% Based on Mätzler, C., & Wegmuller, U. (1987). Dielectric properties of freshwater
% ice at microwave frequencies. *Journal of Physics D: Applied Physics*, 20(12), 1623-1630.
%
% :param frequency: frequency in Hz
% :param temperature: temperature in K
% :raises Exception: if liquid water > 0 or salinity > 0 (model unsuitable)
%     :returns: complex permittivity of pure ice
%     """
FREEZING_POINT = 273.15;
DENSITY_OF_ICE = 916.7;
tempC = temperature - FREEZING_POINT;

if any(tempC < 0)
    error("The water temperature must be higher or equal to {FREEZING_POINT}K")
end

theta = 1 - 300.0 / temperature;

e0 = 77.66 - 103.3 * theta;
e1 = 0.0671 * e0;

f1 = 20.2 + 146.4 * theta + 316 * theta^2;
e2 = 3.52 + 7.52 * theta;
% #  % version of Liebe MPM 1993 uses: e2=3.52
f2 = 39.8 * f1;

Ew = e2 + (e1 - e2) / complex(1, -freqGHz / f2) + (e0 - e1) / complex(1, -freqGHz / f1);

end
