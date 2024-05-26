function [M] = MshockforPratio(Pratio,gam)

% This function calculates the mach number given 
% a pressure ratio across a shock
% [M] = MShockforPratio(P2/P1,gamma)

M = sqrt((Pratio*(gam+1)+gam-1)/(2*gam));

