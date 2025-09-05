function [epsm_r,epsm_i] = TVB_diel_mix(eps_i, eps_h, shape, vi)
%%%Reference:Microwave Radar and Radiometric Remote Sensing, Ulaby, 2014,
%%%University of Michigan Press
%%% This function calculates the dielectric of mixtures of solid. For
%%% simplicity it considered three different shape of inclusions
%%%% 1. Circular Disc Inclusions
%%%% 2. Spaherical Inclusions
%%%% 3. Needle Inclusions
%Input Variables:
    %eps_i: complex dielectric constant of inclusion material
    %eps_h: complex dielectric constant of host material
    %shape: shape of the inclusion
        % 1: circular disc
        % 2: spherical 
        % 3: needle
    % vi: inclusion volume fraction
    
%Outputs:
    %eps_m: complex dielectric constant of mixture
%% Author: Mohammad Mousavi (MM), Dec 2020, NASA-JPL
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  



if shape == 1 %case of thin circular disc inclusions
    eps_m = eps_h + vi/3*(eps_i - eps_h)*(2*eps_i*(1-vi)+eps_h*(1+2*vi))...
        ./(vi*eps_h +(1-vi)*eps_i);

elseif shape ==2 % case of spherical inclusions
    eps_m = eps_h + 3*vi*eps_h*(eps_i -eps_h)./((2*eps_h+eps_i)-vi*(eps_i-eps_h));

elseif shape == 3 % case of needle inclusions
    eps_m = eps_h + vi/3*(eps_i-eps_h)*(eps_h*(5+vi)+(1-vi)*eps_i)./ ...
        (eps_h*(1+ vi)+eps_i*(1-vi));
end

epsm_r=real(eps_m);
epsm_i=abs(imag(eps_m));
end

