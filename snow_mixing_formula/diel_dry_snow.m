function [epsr_ds, epsi_ds] = diel_dry_snow(T, Ps, f )
%%%Reference:Microwave Radar and Radiometric Remote Sensing, Ulaby, 2014,
%%%University of Michigan Press
%Input Variables:
    %T: Temperature in C
    %Ps: Dry Snow Density in g/cm^3
    %f: frequency in GHz
%Output Products:
    %epsr: real part of relative dielectric constant
    %epsi: imaginary part of relative dielectric constant
%author: Mohammad Mousavi (MM), April 2020, NASA-JPL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vi=  Ps /0.9167 ;

[~, epsi_ice] = diel_pure_ice(T,f);

if vi <= 0.45
    epsr_ds = 1 + 1.4667 .* vi + 1.435 .* vi.^3; % 0<vi<0.45
end
if vi > 0.45
    epsr_ds = (1+ 0.4759 .*vi).^3; % vi>0.45
end

epsi_ds = 0.34 * vi * epsi_ice ./(1- 0.42 * vi).^2;


end

