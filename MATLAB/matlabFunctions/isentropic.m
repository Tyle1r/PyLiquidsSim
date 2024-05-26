function [ratios] = isentropic(M,gam)

% This function calculates all the isentropic relations
% call this function like this
% ratios = isentropic(Mach Number, Specific Heat Ratio)
% The return values are
% ratios.P0_P = P0/P = (1 +(gam-1)/2*M^2)^(gam/(gam-1))
% ratios.T0_T = T0/T = 1 +(gam-1)/2*M^2
% ratios.rho0_rho = rho0/rho = (1 +(gam-1)/2*M^2)^(1/(gam-1))
% ratios.A_Astar = A/A^* = 1/M*(2/(1+gam)*(1+(gam-1)/2*M^2)^((gam+1)/(2*(gam-1)))

ratios.P0_P = (1 +(gam-1)/2*M.^2).^(gam/(gam-1));
ratios.T0_T = 1 +(gam-1)/2*M.^2;
ratios.rho0_rho = (1 +(gam-1)/2*M.^2).^(1/(gam-1));
ratios.A_Astar = 1./M.*(2/(1+gam)*(1+(gam-1)/2*M.^2)).^((gam+1)/(2*(gam-1)));

end