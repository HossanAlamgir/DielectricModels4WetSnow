function [epsr, epsi] = diel_brine(T,f,Sb)
%%%Reference:Microwave Radar and Radiometric Remote Sensing, Ulaby, 2014,
%%%University of Michigan Press
%%% This function calculates the dielctric of brine
%Input Variables:
    %T: Temperature in degree C
    %f: frequency in GHz
%Output Products:
    %epsr: real part of relative dielectric constant
    %epsi: imaginary part of relative dielectric constant
%% Author: Mohammad Mousavi (MM), April 2020, NASA-JPL
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
f=f*1e9;
eps0=8.854e-12;
epsInf=4.9;
Delta=25-T;

% Sb = brine_salinity(T);
Nb=Sb.*( (1.707e-2)+(1.205e-5.*Sb)+(4.058e-9.*(Sb.^2)) );

c1=1.0-(1.96e-2.*Delta)+(8.08e-5.*(Delta.^2))-((Nb.*Delta).*(3.02e-5+(3.92e-5.*Delta)+...
    Nb.*(1.72e-5-(6.58e-6.*Delta))));

sigmab_25=Nb.*(10.39-(2.378.*Nb)+(0.683.*(Nb.^2))-(0.135.*(Nb.^3))+...
    (1.01e-2.*(Nb.^4)));

b1=1.0+(0.146e-2.*T.*Nb)-(4.89e-2.*Nb)-(2.97e-2.*(Nb.^2))+(5.64e-3.*(Nb.^3));

taub_0=(1.1109e-10)-(3.824e-12.*T)+(6.938e-14.*(T.^2))-(5.096e-16.*(T.^3));

a1=1.0-(0.255.*Nb)+(5.15e-2.*(Nb.^2))-(6.891e-3.*(Nb.^3));

epsb0_0=88.045-(0.4147.*T)+(6.295e-4.*(T.^2))+(1.075e-5.*(T.^3));


epsb0=epsb0_0.*a1;
taub=taub_0.*b1;
sigmab=sigmab_25.*c1;


eps_b=epsInf + ((epsb0-epsInf)/(1-1i.*f.*taub))+ 1i.*(sigmab/(2.*pi.*f.*eps0));

% epsr=real(eps_b);
% epsi=imag(eps_b);
% 
epsr=epsInf + (epsb0-epsInf)/(1+(f.*taub).^2);
epsi=((f.*taub).*( (epsb0-epsInf)/(1+(f.*taub).^2) )) + ...
    (sigmab/(2.*pi.*f.*eps0));


end