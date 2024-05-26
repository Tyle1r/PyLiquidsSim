function [Ms] = fannoMsforL(fL_D,gam)

% This function calculates the subsonic and supersonic Mach numbers
% associated with a given value of f L / D where f is the friction factor
% L is the length until the flow chokes, and D is the pipe diameter
% Call it like this
% Ms = fannoMsforL(fL_D,gamma)


% Subsonic Root
M=0.02;
for iter = 1:100
    Mfunc = 1+(gam-1)/2*M.^2;
    M1func = 1+(gam-1)/2;
    Eq = M1func/gam*log(M1func./Mfunc)-(1-1./M.^2)/gam +M1func/gam*log(M.^2);
    Err = Eq -fL_D;
    dErrdM = ((gam + 1))/(M*gam) - 2/(M^3*gam) - (M*(gam - 1)*(gam + 1)/2)/(gam*((gam - 1)*M^2/2 + 1));
    M = M - Err/dErrdM;
    M = abs(M);
    if (abs(Err) < 1.0e-5)
        break;
    end
end
if iter > 99 || M > 1
    disp(['did not converge to subsonic fanno flow solution, fL_D = ' num2str(fL_D)]);
    Ms(1) = nan;
else
    Ms(1) = M;
end

% Supersonic Root
M=100;
Mfunc = 1+(gam-1)/2*M.^2;
M1func = 1+(gam-1)/2;
fL_Dmax = M1func/gam*log(M1func./Mfunc)-(1-1./M.^2)/gam +M1func/gam*log(M.^2);
if (fL_D > fL_Dmax)
    % disp('Can not reach that value of fL_D for supersonic case')
    Ms(2) = nan;
    return
end

M= 1+fL_D;
for iter = 1:100
    Mfunc = 1+(gam-1)/2*M.^2;
    M1func = 1+(gam-1)/2;
    Eq = M1func/gam*log(M1func./Mfunc)-(1-1./M.^2)/gam +M1func/gam*log(M.^2);
    Err = Eq -fL_D;
    dErrdM = ((gam + 1))/(M*gam) - 2/(M^3*gam) - (M*(gam - 1)*(gam + 1)/2)/(gam*((gam - 1)*M^2/2 + 1));
    M = M - Err/dErrdM;
    M = abs(M);
    if (abs(Err) < 1.0e-5)
        break;
    end
end
if iter > 99 || M < 1
    % disp(['did not converge to supersonic fanno flow solution, fL_D = ' num2str(fL_D)]);
    Ms(2) = nan;
else
    Ms(2) = M;
end



