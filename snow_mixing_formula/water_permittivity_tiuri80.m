function Ew = water_permittivity_tiuri80(freqGHz, temperature)
%     """Calculates the complex water dielectric constant reported by:
%
%     Tiuri, M. and Schultz, H., Theoretical and experimental studies of microwave radiation from a natural snow field. In Rango, A. , ed.
%     Microwave remote sensing of snowpack properties. Proceedings of a workshop ... Fort Collins, Colorado, May 20-22, 1980.
%     Washington, DC, National Aeronautics and Space Center, 225-234. (Conference Publication 2153.)
%
%     https://ntrs.nasa.gov/api/citations/19810010984/downloads/19810010984.pdf
%
% """


%Constatnts
DENSITY_OF_ICE = 916.7;
DENSITY_OF_WATER = 1000;
FREEZING_POINT = 273.15;

tempC = temperature - FREEZING_POINT;

if any(tempC < 0)
    error("The water temperature must be higher or equal to {FREEZING_POINT}K")
end

e2 = 4.903e-2;

e1 = 87.74 - 0.4008 * tempC + 9.398e-4 * tempC^2 + 1.410e-6 * tempC^3;

% version of Liebe 1991 because Tiuri 1980 does not give the relaxation frequency
theta = 1 - 300.0 / temperature;
f1 = 20.2 + 146.4 * theta + 316 * theta^2;


Ew = e2 + (e1 - e2) / complex(1, -freqGHz /f1);