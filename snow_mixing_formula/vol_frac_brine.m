function vb = vol_frac_brine(T,Si)
%%%Reference:Microwave Radar and Radiometric Remote Sensing, Ulaby, 2014,
%%%University of Michigan Press
%%% This function calculates the volume fraction of brine in sea ice
%Input Variables:
    %T: Temperature in degree C ( -0.5 C >= T >= -22.9 C)
%Output Products:
    %Sb: brine salinity in practical salinity unit (psu)
    
%% Author: Mohammad Mousavi (MM), Dec 2020, NASA-JPL (mousavi@jpl.nasa.gov)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
vb=(1e-3).*Si.*((-49.185./T)+0.532);
end

