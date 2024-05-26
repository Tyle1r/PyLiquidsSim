function [ratios] = shock(M1,gam)

% This function calculates all the shock ratios for
% a given incoming Mach number M1
% call this function like this
% [table] = shock(M1,gam)
% The return values are
% ratios.M2
% ratios.P2_P1 
% ratios.T2_T1
% ratios.rho2_rho1
% ratios.P02_P01
% ratios.Astar2_Astar1

M2 = sqrt((1+(gam-1)/2*M1^2)/(gam*M1^2-(gam-1)/2));
ratios.M2 = M2;
ratios.P2_P1 = (1+gam*M1^2)/(1+gam*M2^2);
ratios.T2_T1 = (1+(gam-1)/2*M1^2)/(1+(gam-1)/2*M2^2);
ratios.rho2_rho1 = table.P2_P1/table.T2_T1;
ratios.P02_P01 = table.P2_P1/table.T2_T1^(gam/(gam-1));
ratios.Astar2_Astar1 = M2/M1*(table.T2_T1)^((gam+1)/(2*(gam-1)));

end