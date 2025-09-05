function [Na,Nb,Nc] = depolarizarion_factor(raspec_ca)
% This function calculates the depolarization factor of a special case
% of spheroid with respect to the axis of rotation, Nc, as a function of
% aspect ratio n=c/a=c/b
% Ref: Jones and Friedman, 2000
% Coder: Alamgir Hossan, JPL, 3/29/2025

Nc = 1/(1+1.6*raspec_ca+0.4*raspec_ca^2);
Na = 0.5*(1-Nc);
Nb = Na;
end