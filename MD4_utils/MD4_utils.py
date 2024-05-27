# Function Adapted From Dr.Brian Helenbrook
import utils.gasDynamics_utils

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
    M_exit, M_super = utils.gasDynamics_utils.MforAratio(Aexit/Astar, gamb)
    isentropic_sub = utils.gasDynamics_utils.isentropic(M_exit,gamb)
    isentropic_super = utils.gasDynamics_utils.isentropic(M_super,gamb)
    # Calculating various Cases
    P_exit_sub_chocked = P0/isentropic_sub['P0_P']
    P_exit_super = P0/isentropic_super['P0_P']

    # CHECK CASE WITH NORMAL SHOCK AT EXIT PLANE OF NOZZLE
    shock = utils.gasDynamics_utils.shock(M_super,gamb)
    P0_postShock = P0*shock['P02_P01']
    M_postShock= shock['M2']
    isentropic_postShock = utils.gasDynamics_utils.isentropic(M_postShock,gamb)

    P_preShock = P_exit_super
    P_postShock = P0_postShock/isentropic_postShock['P0_P']
    if (verbose):
        print('Pressure ratio for just barely choked then subsonic: ' + str(P0 / P_exit_sub_chocked))
        print('Pressure ratio for supersonic matched exit ' +str(P0/P_exit_super))
        print('Pressure ratio for shock at exit ' +str(P0/P_postShock))

    if (Pb > P_exit_sub_chocked):
        if (verbose):
            print('does not choke at nozzle')
        M_exit = utils.gasDynamics_utils.MforPratio(P0/Pb,gamb)
        P_exit = Pb

    elif (Pb < P_exit_super):
        if (verbose):
            print('flow is supersonic exiting nozzle');
            print('Expansion to lower pressure outside of nozzle')
        M_exit = M_super
        P_exit = P0 / isentropic_super['P0_P']
    elif (Pb < P_postShock):
        # Patm is not large enough to keep shock in nozzle
        if (verbose):
            print('flow is supersonic exiting nozzle');
            print('then have oblique shocks after nozzle');
        Mexit = M_super
        Pexit = P0 / isentropic_super['P0_P']
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
            _ , M_preShock = utils.gasDynamics_utils.MforAratio(Ashock/Astar,gamb) # Take supersonic root for flow in front of shock
            # Calculate ratios across shock
            shock = shock(M_preShock,gamb)
            # Calculate stagnation pressure downstream of shock
            P0_postShcok = P0*shock['P02_P01']
            # Calculate new choked throat area downstream of shock
            Astar2 = Astar*shock['Astar2_Astar1']
            # Use exit to throat area ratio to calculate exit mach number
            M_exit , _ = utils.gasDynamics_utils.MforAratio(Aexit/Astar2,gamb) # Take subsonic root because flow downstream is subsonic
            # Calculate ratios at exit
            isentropic_exit = utils.gasDynamics_utils.isentropic(M_exit,gamb)
            # Calculate exit pressure
            P_exit = P0_postShcok/isentropic_exit['P0_P']

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
    isentropic_super = utils.gasDynamics_utils.isentropic(M_exit,gamb)
    T_exit = Tf/isentropic_super['T0_T']
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

    Returns:
    mdot (float): mass flow rate
    thrust (float): thrust produced
    Pexit (float): Exit Pressure
    '''
    # Assume choked flow 
    verbose = 0
    isen_ratios = utils.gasDynamics_utils.isentropic(1,gamb)
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
    _, M_exit = utils.gasDynamics_utils.MforAratio(Ae/Astar,gamb)# Choose supersonic root

    exit_ratios = utils.gasDynamics_utils.isentropic(M_exit,gamb)
    Pe = Pc/exit_ratios['P0_P']
    Te = Tf/exit_ratios['T0_T']
    Ve = M_exit*np.sqrt(gamb*Rb*Te)
    thrust = mdot*Ve+(Pe-Pb)*Ae
    mdot, thrust, Pexit = CvgDvg(Pc, Pb, Tf, Astar, Ae, gamb, Rb, verbose)
    return mdot, thrust, Pexit
def calcM(Isp,g,Mr,mdot,TotalEnergy_kg):
    '''
    Calculates mass needed to deliver the specified energy to a rocket
    using Newton-Rhapson


    Args:
    Isp (float): Specific Impulse
    g (float): Gravity
    mdot (float): mass flow rate
    TotalEngergy_kg (float): Potential + Kinetic Energy per KG
    
    Returns:
    Mp (float): Mass of Propellant
    '''

    Mp = Mr # Initial Guess
    for iter in range(1,11):
        xi = np.log((Mr+Mp)/Mr)
        R = Isp*g**2*(Isp/2*xi**2-(Mr+Mp)/mdot*xi+Mp/mdot)-TotalEnergy_kg
        dRdMp = Isp*g**2*(-1/mdot*xi+1/mdot +(Isp*xi-(Mr+Mp)/mdot)*1/(Mr+Mp))
        Mp = Mp -R/dRdMp

    return Mp

def flightEvents(t,y,alt):
    '''
    Checks events based off of a wanted altitude 

    Args:
    t
    y
    alt

    Returns:
    value (float): current value
    isterminal (boolean): Is it in freefall
    direction (float): 
    '''
    altitude_reached = y[1] - alt
    direction = 0
    isterminal = 1
    return altitude_reached,isterminal,direction


def rocket(t, y, P0, Agrain, rhof, Tf, Rb, gamb, Astar, Ae, a, n, mr, g, v, Cd, Ar):
    
    """
    Finds rates of change of the state variables (position, velocity, and mass) based on the current state and other parameters.

    Args: 
    t (float): time (s)
    y (float): altitude
    P0 (float): stagnation pressure
    Agrain (float): Grain Area
    rhof (float): Fuel Density
    Tf (float): Combustion Temperature
    Rb (float): Fuel Property
    gamb (float): Specific Heat Ratio
    Astar (float): Throat Area
    Ae (float): Exit Area
    a (float): Speed of sound
    n (floar): Ville's Law Constat
    mr (float): mass of rocket
    g (float): gravity
    v (float): velocity
    Cd (float): drag coeff
    Ar (float): cross-sectional area of rocket

    Returns:
    rhs (numpy array): rates of change of the state variables (position, velocity, and mass) based on the current state and other parameters.
    
    """
    x = y[0]
    v = y[1]
    mp = y[2]
    
    # Assuming StandardConditions function is defined
    T, Pb, rho = utils.gasDynamics_utils.StandardConditions(x)
    
    # Reduction factors for Agrain
    reduction_factors = {
        200: 0.62, 190: 0.64, 180: 0.66, 170: 0.68, 160: 0.7, 150: 0.72, 
        140: 0.74, 130: 0.76, 120: 0.78, 110: 0.8, 100: 0.82, 90: 0.84, 
        80: 0.86, 70: 0.88, 60: 0.9, 50: 0.92, 40: 0.94, 30: 0.96, 20: 0.98
    }
    
    # Apply reduction factor if t exceeds the specified values
    for threshold, factor in reduction_factors.items():
        if t > threshold:
            Agrain *= factor
    
    if mp > 0:
        # Assuming calc_rocket function is defined
        thrust, mdot, rdot = calc_rocket(Pb, Agrain, rhof, Tf, Rb, gamb, Astar, Ae, a, n)
    else:
        thrust = 0
        mdot = 0
    
    if mdot == 0:
        Isp = 0
    else:
        Isp = thrust / (mdot * g)
    

    # Calculate right-hand side of ODEs
    rhs = np.zeros(3)

    # dx/dt = v
    rhs[0] = v

    # dv/dt = F/m
    W = (mr + mp) * g
    D = 0.5 * rho * v**2 * Cd * Ar
    rhs[1] = (thrust - W - D) / (mr + mp)

    # dmp / dt = -mdot
    rhs[2] = -mdot

    return rhs
'''
def rocket_rhs(t, y, P0, Agrain, rhof, Tf, Rb, gamb, Astar, Ae, a, n, mr, g, v, Cd, Ar):
    x = y[0]
    v = y[1]
    mp = y[2]
    
    # Assuming StandardConditions function is defined
    T, Pb, rho = utils.gasDynamics_utils.StandardConditions(x)
    
    # Reduction factors for Agrain
    reduction_factors = {
        200: 0.62, 190: 0.64, 180: 0.66, 170: 0.68, 160: 0.7, 150: 0.72, 
        140: 0.74, 130: 0.76, 120: 0.78, 110: 0.8, 100: 0.82, 90: 0.84, 
        80: 0.86, 70: 0.88, 60: 0.9, 50: 0.92, 40: 0.94, 30: 0.96, 20: 0.98
    }
    
    # Apply reduction factor if t exceeds the specified values
    for threshold, factor in reduction_factors.items():
        if t > threshold:
            Agrain *= factor
    
    if mp > 0:
        # Assuming calc_rocket function is defined
        thrust, mdot, rdot = calc_rocket(Pb, Agrain, rhof, Tf, Rb, gamb, Astar, Ae, a, n)
    else:
        thrust = 0
        mdot = 0

    Isp = thrust / (mdot * g)

    # Calculate right-hand side of ODEs
    rhs = np.zeros(3)

    # dx/dt = v
    rhs[0] = v

    # dv/dt = F/m
    W = (mr + mp) * g
    D = 0.5 * rho * v**2 * Cd * Ar
    rhs[1] = (thrust - W - D) / (mr + mp)

    # dmp / dt = -mdot
    rhs[2] = -mdot

    return rhs
    '''