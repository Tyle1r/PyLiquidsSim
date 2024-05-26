function [Ms] = rayleighsMsforT0(T0_Tstar,gam)

% This function calculates the subsonic & supersonic Mach numbers that result in a given
% ratio of T0/T0* for Rayleigh flow
% Call it like this
% [Ms] = rayleighsMsforT0(T0_Tstar,gam)

% Answer was found using 
% syms M gam T0_Tstar positive real
% solve(T0_Tstar == (1+gam)^2*M.^2./(1+gam*M.^2).^2.*(2+(gam-1)*M.^2)/(gam+1),M)
Ms(1) = (-(T0_Tstar*gam - gam + gam*(1 - T0_Tstar)^(1/2) + (1 - T0_Tstar)^(1/2) - 1)/(T0_Tstar*gam^2 - gam^2 + 1))^(1/2);
Ms(2) = ((gam - T0_Tstar*gam + gam*(1 - T0_Tstar)^(1/2) + (1 - T0_Tstar)^(1/2) + 1)/(T0_Tstar*gam^2 - gam^2 + 1))^(1/2);
