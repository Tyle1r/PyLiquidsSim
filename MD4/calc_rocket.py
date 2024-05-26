# Function Adapted From Dr.Brian Helenbrook
from utils.gasDynamics_utils import isentropic, MforAratio, shock, MforPratio
import numpy as np
def CvgDvg(P0,Pb,Tf,Astar,Aexit,gamb,R,verbose):
    '''
    Calculates the thrust of converging diverging nozzle

    Args:
    P0 (float): Stagnation Pressure [Pa]
    Pb (float): Back Pressure [Pa]
    T0 (float): Stagnation Temperature [K]
    Astar (float): Nozzle Throat Area [m^2]
    gam (float): Specific Heat Ratio
    R (float): Gas Constant
    verbose (boolean): Output Information About Nozzle Conditions
    '''
    verbose = 0
    Mstar = 1
    M_exit, M_super = MforAratio(Aexit/Astar, gamb)
    isentropic_sub = isentropic(M_exit,gamb)
    isentropic_super = isentropic(M_super,gamb)
    # Calculating various Cases
    P_exit_sub_chocked = P0/isentropic_sub.P0_P
    P_exit_super = P0/isentropic_super.P0_P

    # CHECK CASE WITH NORMAL SHOCK AT EXIT PLANE OF NOZZLE
    shock = shock(M_super,gamb)
    P0_postShock = P0*shock.P02_P01
    M_postShock= shock.M2
    isentropic_postShock = isentropic(M_postShock,gamb)

    P_preShock = P_exit_super
    P_postShock = P0_postShock/isentropic_postShock.P0_P
    if (verbose):
        print('Pressure ratio for just barely choked then subsonic: ' + str(P0 / P_exit_sub_chocked))
        print('Pressure ratio for supersonic matched exit ' +str(P0/P_exit_super))
        print('Pressure ratio for shock at exit ' +str(P0/P_postShock))

    if (Pb > P_exit_sub_chocked):
        if (verbose):
            print('does not choke at nozzle')
        M_exit = MforPratio(P0/Pb,gamb)
        P_exit = Pb

    elif (Pb < P_exit_super):
        if (verbose):
            print('flow is supersonic exiting nozzle');
            print('Expansion to lower pressure outside of nozzle')
        M_exit = M_super
        P_exit = P0 / isentropic_super.P0_P
    elif (Pb < P_postShock):
        # Patm is not large enough to keep shock in nozzle
        if (verbose):
            print('flow is supersonic exiting nozzle');
            print('then have oblique shocks after nozzle');
        Mexit = M_super
        Pexit = P0 / isentropic_super.P0_P
    else:
        # Use bisection method to determine shock location
        # For book problems where shock location is given, only
        # need to do this once (no iterating).
        # Guess shock is in the middle of my two values
        A_before = Astar
        A_after = Aexit

    for iter in range(1,101):
        Ashock = 0.5*(A_before+A_after)

        # Find mach number right before shock
        _ , M_preShock = MforAratio(Ashock/Astar,gamb) # Take supersonic root for flow in front of shock
        # Calculate ratios across shock
        shock = shock(M_preShock,gamb)
        # Calculate stagnation pressure downstream of shock
        P0_postShcok = P0*shock.P02_P01
        # Calculate new choked throat area downstream of shock
        Astar2 = Astar*shock.Astar2_Astar1
        # Use exit to throat area ratio to calculate exit mach number
        M_exit , _ = MforAratio(Aexit/Astar2,gamb) # Take subsonic root because flow downstream is subsonic
        # Calculate ratios at exit
        isentropic_exit = isentropic(M_exit,gamb)
        # Calculate exit pressure
        P_exit = P0_postShcok/isentropic_exit.P0_P

        if (Pexit > Pb):
            # Exit pressure is too high
            # Didn't expand enough
            # Shock is too far up throat
            Aleft = Ashock
        else:
            # Shock is too far down throat
            Aright = Ashock
    if (verbose):
        print('Shock is inside nozzle')
        print('Shock area ratio is: ' +str(Ashock/Astar))
        print('Shock Mach number is: ' +str(M_preShock))
    isentropic_super = isentropic(M_exit,gamb)
    T_exit = Tf/isentropic_super.T0_T
    c_exit = np.sqrt(gamb*R*T_exit)
    v_exit = c_exit*M_exit
    rho_exit = Pexit/(R*T_exit)
    m_dot = rho_exit*v_exit*Aexit
    thrust = m_dot*v_exit +(Pexit-Pb)*Aexit
    if (verbose):
        print('Pressure Ratio P0/Pexit: ' +str(P0/Pexit))
        print('Exit Mach Number: ' +str(M_exit))
        print('Exit Temperature: ' +str(T_exit))
        print('Exit velocity: ' +str(v_exit))
    return m_dot,thrust,Pexit

def calc_rocket(Pb,Agrain,rhof,Tf,Rb,gamb,Astar,Ae,a,n):
    '''
    This Function calculates the thrust, mass flow, and rate of combustion
    for a solid rocket motor

    Args:
    Pb (float): Back Pressure (Pa)
    Agrain (float): Solid Grain X-Area (m^2)
    rhof (float): fuel density (kg/m^3)
    Tf (float): Combustion Temperature (K)
    Rb (float): Fuel Thermo Constant
    gamb (float): Specific Heat Ratio
    Astar (float): Throat of Nozzle Area (m^2)
    Ae (float): Exit Nozzle Area (m^2)
    a (float): Speed of Sound (m/s)
    n (float): Vielle’s law Burn Constant
    '''
    # Assume choked flow 
    verbose = 0
    isen_ratios = isentropic(1,gamb)
    Tstar = Tf/((gamb+1)/2)
    Vstar = np.sqrt(gamb*Rb*Tstar)
    '''
    Find chamber stagnation pressure
    Need burning rate mass flow rate to be equal to 
    flow rate through throat
    
    Use rdot*Agrain*rhof = rho_star*Vstar*Athroat
    to the stagnation pressure in the rocket
    '''
    Pc = (1/((gamb+1)/2)**(gamb/(gamb-1))*1/(Rb*Tstar)*Vstar/(rhof*a)*Astar/Agrain)**(1/(n-1))
    # Calculate mdot based on burn rate (Vielle’s law)
    rdot = a*Pc**n
    mdot = rhof*rdot*Agrain

    # Calculate mdot based on throat conditions
    # Double checking that these come out equal

    Pstar = Pc/((gamb+1)/2)**(gamb/(gamb-1))
    rhostar = Pstar/(Rb*Tstar)
    mdot = Pc/((gamb+1)/2)**(gamb/(gamb-1))*1/(Rb*Tstar)*Astar*Vstar

    # Now calculate exit conditions
    Ms = MforAratio(Ae/Astar,gamb)
    Me = MforAratio.M_super; # Choose supersonic root

    exit_ratios = isentropic(Me,gamb)
    Pe = Pc/exit_ratios.P0_P
    Te = Tf/exit_ratios.T0_T
    Ve = Me*np.sqrt(gamb*Rb*Te)
    thrust = mdot*Ve+(Pe-Pb)*Ae
    mdot, thrust, Pexit = CvgDvg(Pc, Pb, Tf, Astar, Ae, gamb, Rb, verbose)