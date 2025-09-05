function Sb = brine_salinity(T)
%%%Reference:Microwave Radar and Radiometric Remote Sensing, Ulaby, 2014,
%%%University of Michigan Press
%Input Variables:
    %T: Temperature in degree C
%Output Products:
    %Sb: brine salinity in practical salinity unit (psu)
    
%% Author: Mohammad Mousavi (MM), Dec 2020, NASA-JPL (mousavi@jpl.nasa.gov)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
if T>= -8.2 && T<=-2
    Sb=1.725 - (18.756.*T) - (0.3964.*(T.^2));
elseif T >= -22.9 && T<= -8.2
    Sb=57.041 - (9.929.*T) - (0.16204.*(T.^2)) - (0.002396.*(T.^3));
elseif T>= -36.8 && T<= -22.9
    Sb=242.94 + (1.5299.*T) + (0.0429.*(T.^2));
elseif T>= -43.2 && T<= -36.8 
    Sb=508.18 + (14.535.*T) + (0.2018.*(T.^2));
end
end

