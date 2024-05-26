function [M] = MforPratio(Pratio,gam)

% This function calculates the mach number given a stagnation to static pressure ratio
% [M] = MforPratio(P0/P,gamma)

M = sqrt((Pratio^((gam-1)/gam)-1)*2/(gam-1));