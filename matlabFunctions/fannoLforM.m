function [data] = fannoLforM(M,gam)

% This function calculates the the value of f L / D as well as ratios to
% sonic conditions associated with a given value of M where f is the friction factor
% L is the length until the flow chokes, and D is the pipe diameter
% Call it like this
% data = fannoLforM(M,gamma)

Mfunc = 1+(gam-1)/2*M.^2;
M1func = 1+(gam-1)/2;

data.fL_D = M1func/gam*log(M1func./Mfunc)-(1-1./M.^2)/gam +M1func/gam*log(M.^2);
data.T1_Tstar = M1func./Mfunc;
data.P1_Pstar = sqrt(M1func./Mfunc)./M;
data.P0_P0star = data.P1_Pstar.*(Mfunc./M1func).^(gam/(gam-1));
data.rho_rhostar = data.P1_Pstar./data.T1_Tstar;

