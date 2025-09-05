function [epsr, epsi] = diel_pure_ice(T,f)
%%%Reference:Microwave Radar and Radiometric Remote Sensing, Ulaby, 2014,
%%%University of Michigan Press
%Input Variables:
    %T: Temperature in degree C
    %f: frequency in GHz
%Output Products:
    %epsr: real part of relative dielectric constant
    %epsi: imaginary part of relative dielectric constant
%author: Mohammad Mousavi (MM), April 2020, NASA-JPL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

T = T + 273; 

theta = (300./T) - 1;

B1 = 0.0207;
B2 = 1.16e-11;
b = 335;

alpha = (0.00504 + 0.0062.*theta).*exp(-22.1*theta);

betaM = (B1./T) .* exp(b./T) ./ (exp(b./T)-1).^2  + B2.*f.^2;
delBeta = exp(-9.963 + 0.0372.*(T-273.16));

beta = betaM + delBeta;

epsr = 3.1884 + 9.1e-4 *(T-273);

epsi = alpha./f + beta.*f;
end