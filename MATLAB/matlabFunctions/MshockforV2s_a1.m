function [M] = MshockforV2s_a1(V2s_a1,gam)

% This function calculates the mach number given 
% the ratio of the downstream velocity to the upstream speed of sound
% It is assumed that fluid upstream of the shock is stationary
% [M] = MforV2s(V2s/a1,gamma)

M = (gam+1)/4*abs(V2s_a1)+sqrt(((gam+1)/4*V2s_a1)^2+1);