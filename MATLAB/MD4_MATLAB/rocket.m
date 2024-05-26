function [rhs] = rocket(t,y,P0,Agrain,rhof,Tf,Rb,gamb,Astar,Ae,a,n,mr,g,v,Cd,Ar)

x = y(1);
v = y(2);
mp = y(3);
[T,Pb,rho] = StandardConditions(x);
%Keep Agrain constant until Gs exceed


% if t > 260
%     Agrain = 0.5 * Agrain; % Reduce Agrain by 10% in the first 10 seconds
% elseif t > 250
%     Agrain = 0.52 * Agrain; % Reduce Agrain by 20% from 10 to 20 seconds
% elseif t > 240
%     Agrain = 0.54 * Agrain; % Reduce Agrain by 20% from 10 to 20 seconds
% elseif t > 230
%     Agrain = 0.56 * Agrain; % Reduce Agrain by 20% from 10 to 20 seconds
% elseif t > 220
%     Agrain = 0.58 * Agrain; % Reduce Agrain by 20% from 10 to 20 seconds
% elseif t > 210
%     Agrain = 0.6 * Agrain; % Reduce Agrain by 20% from 10 to 20 seconds
if t > 200
    Agrain = 0.62 * Agrain; % Reduce Agrain by 20% from 10 to 20 seconds
elseif t > 190
    Agrain = 0.64 * Agrain; % Reduce Agrain by 20% from 10 to 20 seconds
elseif t > 180
    Agrain = 0.66 * Agrain; % Reduce Agrain by 20% from 10 to 20 seconds
elseif t > 170
    Agrain = 0.68 * Agrain; % Reduce Agrain by 20% from 10 to 20 seconds
elseif t > 160
    Agrain = 0.7 * Agrain; % Reduce Agrain by 20% from 10 to 20 seconds
elseif t > 150
    Agrain = 0.72 * Agrain; % Reduce Agrain by 20% from 10 to 20 seconds
elseif t > 140
    Agrain = 0.74 * Agrain; % Reduce Agrain by 20% from 10 to 20 seconds
elseif t > 130
    Agrain = 0.76 * Agrain; % Reduce Agrain by 20% from 10 to 20 seconds
elseif t > 120
    Agrain = 0.78 * Agrain; % Reduce Agrain by 20% from 10 to 20 seconds
elseif t > 110
    Agrain = 0.8 * Agrain; % Reduce Agrain by 20% from 10 to 20 seconds
elseif t > 100
    Agrain = 0.82 * Agrain; % Reduce Agrain by 20% from 10 to 20 seconds
elseif t > 90
    Agrain = 0.84 * Agrain; % Reduce Agrain by 20% from 10 to 20 seconds
elseif t > 80
    Agrain = 0.86 * Agrain; % Reduce Agrain by 20% from 10 to 20 seconds
elseif t > 70
    Agrain = 0.88 * Agrain; % Reduce Agrain by 20% from 10 to 20 seconds
elseif t > 60
    Agrain = 0.9 * Agrain; % Reduce Agrain by 20% from 10 to 20 seconds
elseif t > 50
    Agrain = 0.92 * Agrain; % Reduce Agrain by 20% from 10 to 20 seconds
elseif t > 40
    Agrain = 0.94 * Agrain; % Reduce Agrain by 20% from 10 to 20 seconds
elseif t > 30
    Agrain = 0.96 * Agrain; % Reduce Agrain by 20% from 10 to 20 seconds
elseif t > 20
    Agrain = 0.98 * Agrain; % Reduce Agrain by 20% from 10 to 20 seconds
end

if mp > 0
    [thrust,mdot,rdot] = calc_rocket(Pb,Agrain,rhof,Tf,Rb,gamb,Astar,Ae,a,n);
else
    thrust = 0;
    mdot = 0;
end

Isp = thrust/(mdot*g);

rhs = zeros(3,1);

% dx/dt = v
rhs(1) = v;

% dv/dt = F/m
W = (mr+mp)*g;
D = 1/2*rho*v^2*Cd*Ar;
rhs(2) = (thrust-W-D)/(mr+mp);

% dmp / dt = mdot
rhs(3) = -mdot;

end
