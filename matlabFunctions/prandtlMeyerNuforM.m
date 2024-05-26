function [nu] = prandtlMeyerNuforM(M,gam)

% Script to calculate Prandtl-Meyer Function
% inputs are value of M and gamma
% output is the Prantdtl-Meyer function, nu
% Call it like this
% [nu] = prandtlMeyerNuforM(M,gam)

const = sqrt((gam+1)/(gam-1));
nu = const*atan(sqrt(M^2-1)/const) -atan(sqrt(M^2-1));
    

