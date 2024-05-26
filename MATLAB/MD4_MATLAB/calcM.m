function [Mp] = calcM(Isp,g,Mr,mdot,TotalEnergy_kg)
% [Mass of Propellant] = calcM(Isp,g,Mr,mdot,TotalEnergy_kg)
% The function uses Newton-Rhapson to calculate the mass needed
% To deliver the specified energy to the rocket
% Based on the following equation

%% 
% $E/kg = I_{sp} g^2 \left( \left(I_{sp} \frac{1}{2} \log\left(\frac{M_r+M_p}{M_r}\right) 
% -\frac{M_r+M_p}{\dot{m}} \right) \log\left(\frac{M_r+M_p}{M_r}\right) +\frac{M_p}{\dot{m}}\right)$ 
% to estimate 


% Isp is the specific impulse of the rocket
% g is the acceleration of gravity
% Mr is the mass of the rocket (excluding propellant)
% mdot is the mass flow rate out of the rocket
% TotalEnergy_kg is the desired potential+kinetic energy per kg 
% to be delivered to the rocket after burning is complete

% Function return Mp which is the mass of the propellant

% Initial Guess
Mp = Mr;
for iter = 1:10
        xi = log((Mr+Mp)/Mr);
        R = Isp*g^2.*(Isp/2.*xi.^2-(Mr+Mp)./mdot.*xi+Mp./mdot)-TotalEnergy_kg;
        dRdMp = Isp*g^2.*(-1./mdot.*xi+1./mdot +(Isp.*xi-(Mr+Mp)./mdot).*1./(Mr+Mp));
        Mp = Mp -R./dRdMp;
    end
end