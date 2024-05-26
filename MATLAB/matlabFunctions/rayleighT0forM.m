function [data] = rayleighT0forM(M,gam)
% This calculates ratios of P/P*, T/T*, T0/T* and V/V*
% For Rayleigh flow for a given Mach number
% Call it like this
% [data] = rayleighT0forM(M,gam)


data.P_Pstar = (1+gam)./(1+gam*M.^2);
data.T_Tstar = (1+gam)^2*M.^2./(1+gam*M.^2).^2;
data.T0_T0star = data.T_Tstar.*(2+(gam-1)*M.^2)/(gam+1);
data.P0_P0star = data.P_Pstar.*((2+(gam-1)*M.^2)/(gam+1)).^(gam/(gam-1));
data.V_Vstar = sqrt(data.T_Tstar).*M;