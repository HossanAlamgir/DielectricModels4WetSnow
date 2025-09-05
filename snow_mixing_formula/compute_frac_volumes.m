function [frac_volume, fi, fw] =  compute_frac_volumes(density, liquid_water)
%     """Computes the fractional volume of ice+water, the fractional volume of ice, and the fractional volume of water
%     from the (wet) snow density and the liquid_water which is the volume fraction of liquid with respect to ice + liquid (but no air).
%
%     :param density: density of the snow, including the ice and water phases.
%     :param liquid_water: (fractional volume of water with respect to ice+water volume).
%
%     :returns: frac_volume, fi, fw
% """

%Constatnts
DENSITY_OF_ICE = 916.7;
DENSITY_OF_WATER = 1000;
FREEZING_POINT = 273.15;
density_melange = DENSITY_OF_ICE * (1 - liquid_water) + DENSITY_OF_WATER * liquid_water;
% # variations of density with temperature (a few kg/m3) and air mass (less than 1 kg/m3) are not taken into account.

frac_volume = density / density_melange;

fi = frac_volume * (1 - liquid_water);

fw = frac_volume * liquid_water;

end