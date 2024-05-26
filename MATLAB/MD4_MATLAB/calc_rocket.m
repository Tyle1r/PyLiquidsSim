function [thrust,mdot,rdot] = calc_rocket(Pb,Agrain,rhof,Tf,Rb,gamb,Astar,Ae,a,n)
% Assume choked flow 
verbose = 0;
ratios = isentropic(1,gamb);
Tstar = Tf/((gamb+1)/2);
Vstar = sqrt(gamb*Rb*Tstar);

% Find chamber stagnation pressure
% Need burning rate mass flow rate to be equal to 
% flow rate through throat

% Use rdot*Agrain*rhof = rho_star*Vstar*Athroat
% to the stagnation pressure in the rocket
Pc = (1/((gamb+1)/2)^(gamb/(gamb-1))*1/(Rb*Tstar)*Vstar/(rhof*a)*Astar/Agrain)^(1/(n-1));

% Calculate mdot based on burn rate
rdot = a*Pc^n;
mdot = rhof*rdot*Agrain;

% Calculate mdot based on throat conditions
% Double checking that these come out equal
Pstar = Pc/((gamb+1)/2)^(gamb/(gamb-1));
rhostar = Pstar/(Rb*Tstar);
mdot = Pc/((gamb+1)/2)^(gamb/(gamb-1))*1/(Rb*Tstar)*Astar*Vstar;

% Now calculate exit conditions
Ms = MforAratio(Ae/Astar,gamb);
Me = Ms(2); % Choose supersonic root

exit_ratios = isentropic(Me,gamb);
Pe = Pc/exit_ratios.P0_P;
Te = Tf/exit_ratios.T0_T;
Ve = Me*sqrt(gamb*Rb*Te);
thrust = mdot*Ve+(Pe-Pb)*Ae;

[mdot,thrust,Pexit] = CvgDvg(Pc,Pb,Tf,Astar,Ae,gamb,Rb,verbose);
