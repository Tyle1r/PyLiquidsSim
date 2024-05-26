function [mdot,thrust,Pexit] = CvgDvg(P0,Pb,Tf,Astar,Aexit,gamb,R,verbose)

% This script calculates the thrust of converging diverging nozzle
% The nozzle is assumed to have 100% efficiency
% but there may be shocks in the nozzle
% The inputs are
% P0 stagnation pressure
% Pb back pressure
% T0 stagnation temperature
% Astar throat area
% Aexit exit area
% gam specific heat ratio
% R gas constant
% verbose (T/F) output information about nozzle conditions
verbose = 0;
Mstar = 1;
Ms = MforAratio(Aexit/Astar, gamb);
Mexit = Ms(1);
ratiossub = isentropic(Mexit,gamb);
Msup = Ms(2);
ratiossuper = isentropic(Msup,gamb);
%Calculating various Cases

pexitsubsonic_choked = P0/ratiossub.P0_P;
pexitsupersonic = P0 / ratiossuper.P0_P;

% CHECK CASE WITH NORMAL SHOCK AT EXIT PLANE OF NOZZLE
data = shock(Msup,gamb);
P02 = P0*data.P02_P01;
M2= data.M2;
ratiosafter = isentropic(M2,gamb);

Pbeforeshock = pexitsupersonic;
Paftershock = P02/ratiosafter.P0_P;
if (verbose)
    disp(['Pressure ratio for just barely choked then subsonic ' num2str(P0/pexitsubsonic_choked)]);
    disp(['Pressure ratio for supersonic matched exit ' num2str(P0/pexitsupersonic)]);
    disp(['Pressure ratio for shock at exit ' num2str(P0/Paftershock)]);
end

if (Pb > pexitsubsonic_choked)
    if (verbose)
        disp('does not choke at nozzle');
    end
    Mexit = MforPratio(P0/Pb,gamb);
    Pexit = Pb;

elseif (Pb < pexitsupersonic)
    if (verbose)
        disp('flow is supersonic exiting nozzle');
        disp('Expansion to lower pressure outside of nozzle')
    end
    Mexit = Msup;
    Pexit = P0 / ratiossuper.P0_P;
elseif (Pb < Paftershock)
    % Patm is not large enough to keep shock in nozzle
    if (verbose)
        disp('flow is supersonic exiting nozzle');
        disp('then have oblique shocks after nozzle');
    end
    Mexit = Msup;
    Pexit = P0 / ratiossuper.P0_P;
else
    % Use bisection method to determine shock location
    % For book problems where shock location is given, only
    % need to do this once (no iterating).
    % Guess shock is in the middle of my two values
    Aleft = Astar;
    Aright = Aexit;

    for iterations=1:100
        Ashock = 0.5*(Aleft+Aright);

        % Find mach number right before shock
        Ms = MforAratio(Ashock/Astar,gamb);
        % Take supersonic root for flow in front of shock
        Mbeforeshock = Ms(2);
        % Calculate ratios across shock
        shock_data = shock(Mbeforeshock,gamb);
        % Calculate stagnation pressure downstream of shock
        P02 = P0*shock_data.P02_P01;
        % Calculate new choked throat area downstream of shock
        Astar2 = Astar*shock_data.Astar2_Astar1;
        % Use exit to throat area ratio to calculate exit mach number
        Ms = MforAratio(Aexit/Astar2,gamb);
        % Take subsonic root because flow downstream is subsonic
        Mexit = Ms(1); 
        % Calculate ratios at exit
        exit_ratios = isentropic(Mexit,gamb);
        % Calculate exit pressure
        Pexit = P02/exit_ratios.P0_P;

        if (Pexit > Pb)
            % Exit pressure is too high
            % Didn't expand enough
            % Shock is too far up throat
            Aleft = Ashock;
        else
            % Shock is too far down throat
            Aright = Ashock;
        end
    end
    if (verbose)
        disp('Shock is inside nozzle');
        disp(['Shock area ratio is: ' num2str(Ashock/Astar)])
        disp(['Shock Mach number is: ' num2str(Mbeforeshock)])
    end
end

ratiossuper = isentropic(Mexit,gamb);
Texit = Tf/ratiossuper.T0_T;
cexit = sqrt(gamb*R*Texit);
vexit = cexit*Mexit;
rhoexit = Pexit/(R*Texit);
mdot = rhoexit*vexit*Aexit;
thrust = mdot*vexit +(Pexit-Pb)*Aexit;

if (verbose)
    disp(['Pressure Ratio P0/Pexit: ' num2str(P0/Pexit)])
    disp(['Exit Mach Number: ' num2str(Mexit)])
    disp(['Exit Temperature: ' num2str(Texit)])
    disp(['Exit velocity: ' num2str(vexit)])
end
end