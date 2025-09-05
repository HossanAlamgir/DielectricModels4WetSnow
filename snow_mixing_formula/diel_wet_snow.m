function [epsr, epsi] = diel_wet_snow(Ps, mv, f )
%%%Reference:Microwave Radar and Radiometric Remote Sensing, Ulaby, 2014,
%%%University of Michigan Press
%Input Variables:
    %Ps: Dry Snow Density (g/cm^3)
    %mv: volumetric water content (0<mv<30 in percentage scale)
    %f: frequency in GHz
%Output Products:
    %epsr: real part of relative permitivity
    %epsi: imaginary part of relative permitivity
%author: Mohammad Mousavi (MM), April 2020, NASA-JPL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A1 =  0.78 + 0.03 .*f - 0.58e-3 .* f.^2;
A2 =  0.97 - 0.39e-2 .*f + 0.39e-3 .* f.^2;
B1 =  0.31 - 0.05 .*f + 0.87e-3 .* f.^2;

A = A1 .*(1.0 + 1.83*Ps + 0.02*mv.^1.015) + B1;
B = 0.073 .*A1;
C = 0.073 .*A2;
x = 1.31;
f0 = 9.07;

epsr = A + (B .* mv.^x) ./ (1+(f/f0).^2);
epsi = (C .* (f/f0) .* mv.^x) ./ (1+(f/f0).^2);
end
