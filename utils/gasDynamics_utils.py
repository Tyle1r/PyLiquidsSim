# All functions are adapted from MATLAB functions written by Dr.Brian Helenbrook

import numpy as np
from scipy.special import cot
def isentropic(M, gam):
    """
    Calculate isentropic flow property ratios for given Mach numbers and specific heat ratio.

    Args:
    M (numpy.ndarray): Array of Mach numbers.
    gam (float): Specific heat ratio.

    Returns:
    isen_ratios: Dictionary containing the ratios P0_P, T0_T, rho0_rho, and A_Astar.
    """
    isen_ratios = {
        'P0_P': (1 + (gam - 1) / 2 * M**2)**(gam / (gam - 1)),
        'T0_T': 1 + (gam - 1) / 2 * M**2,
        'rho0_rho': (1 + (gam - 1) / 2 * M**2)**(1 / (gam - 1)),
        'A_Astar': 1 / M * (2 / (1 + gam) * (1 + (gam - 1) / 2 * M**2))**((gam + 1) / (2 * (gam - 1)))
    }
    return isen_ratios

def fannoLforM(M,gam):
    '''
    Calculate fanno L/D and ratios for sonic conditions
    where f is the friction factor L is the length until the flow chokes, and D is the pipe diameter
    
    Args:
    M (numpy.ndarray): Array of Mach numbers.
    gam (float): Specific heat ratio.
    
    Returns:
    fanno_ratios: Dictionary containing the ratios fL_D, T1_Tstar, P1_Pstar, P0_P0star,rho_rhostar
    '''
    Mfunc = 1 + (gam - 1) / 2 * M**2
    M1func = 1 + (gam - 1) / 2

    fanno_ratios = {
        'fL_D': M1func/gam*np.log(M1func/Mfunc) - (1 - 1/M**2)/gam + M1func/gam*np.log(M**2),
        'T1_Tstar': M1func/Mfunc,
        'P1_Pstar': np.sqrt(M1func/Mfunc) / M
    }
    
    # Adding P0_P0star after P1_Pstar is calculated
    fanno_ratios['P0_P0star'] = fanno_ratios['P1_Pstar'] * (Mfunc / M1func)**(gam / (gam - 1))
    # Adding rho_rhostar after P1_Pstar and T1_Tstar are calculated
    fanno_ratios['rho_rhostar'] = fanno_ratios['P1_Pstar'] / fanno_ratios['T1_Tstar']

    return fanno_ratios


def fannoMsforL(fL_D,gam):
    '''
    This function calculates the subsonic and supersonic Mach numbers
    associated with a given value of f L / D where f is the friction factor
    L is the length until the flow chokes, and D is the pipe diameter

    Args:
    fL_D (numpy.ndarray): Array of Mach numbers.
    gam (float): Specific heat ratio.
    
    Returns:
    Ms_subsonic, Ms_supersonic (float): Subsonic and supersonic Mach numbers
    '''
    # Subsonic root
    M = 0.02
    for iter in range(1, 101):
        Mfunc = 1 + (gam - 1) / 2 * M**2
        M1func = 1 + (gam - 1) / 2
        Eq = M1func / gam * np.log(M1func / Mfunc) - (1 - 1 / M**2) / gam + M1func / gam * np.log(M**2)
        Err = Eq - fL_D
        dErrdM = (gam + 1) / (M * gam) - 2 / (M**3 * gam) - (M * (gam - 1) * (gam + 1) / 2) / (gam * ((gam - 1) * M**2 / 2 + 1))
        M = M - Err / dErrdM
        M = abs(M)
        if abs(Err) < 1.0e-5:
            break

    if iter > 99 or M > 1:
        print(f'did not converge to subsonic fanno flow solution, fL_D = {fL_D}')
        Ms_subsonic = np.NaN
    else:
        Ms_subsonic = M

    # Supersonic root
    M = 100
    Mfunc = 1 + (gam - 1) / 2 * M**2
    M1func = 1 + (gam - 1) / 2
    fL_Dmax = M1func / gam * np.log(M1func / Mfunc) - (1 - 1 / M**2) / gam + M1func / gam * np.log(M**2)
    if fL_D > fL_Dmax:
        # print('Can not reach that value of fL_D for supersonic case')
        Ms_supersonic = np.NaN
        return Ms_subsonic, Ms_supersonic

    M = 1 + fL_D
    for iter in range(1, 101):
        Mfunc = 1 + (gam - 1) / 2 * M**2
        M1func = 1 + (gam - 1) / 2
        Eq = M1func / gam * np.log(M1func / Mfunc) - (1 - 1 / M**2) / gam + M1func / gam * np.log(M**2)
        Err = Eq - fL_D
        dErrdM = (gam + 1) / (M * gam) - 2 / (M**3 * gam) - (M * (gam - 1) * (gam + 1) / 2) / (gam * ((gam - 1) * M**2 / 2 + 1))
        M = M - Err / dErrdM
        M = abs(M)
        if abs(Err) < 1.0e-5:
            break

    if iter > 99 or M < 1:
        # print(f'did not converge to supersonic fanno flow solution, fL_D = {fL_D}')
        Ms_supersonic = np.NaN
    else:
        Ms_supersonic = M

    return Ms_subsonic, Ms_supersonic


def ffactor(Re,eps_D):

    '''
    Calculates pipe friction factor

    Args:
    Re (float): Reynolds number of the pipe flow
    eps_D (float): roughness height divided by pipe diameter

    Returns:
    f (float): Pipe friction factor
    '''


    #Just do a fixed point iteration
    f = 0.02
    for i in range(1, 101):
        fnew = 1/(-2*np.log10(eps_D/3.7 +2.51/(Re*np.sqrt(f))))**2
        if (abs(f-fnew) < 1.0e-8):
            break
        f = fnew
    if i > 99:
        print('Trouble converging to friction factor')
        return f

def MforAratio(Aratio,gam):
    '''
    Calculates the mach number of flow with a given area of expansion

    Args:
    Aratio (float): Area ratio of expansion
    gam (float): Specific Heat Ratio

    Returns:
    M_sub (float): Subsonic Mach
    M_super (float): Supersonic Mach 
    '''
    if Aratio<1:
        print('Area ratio is less than 1: ' + str(Aratio))
        M_sub = 1
        M_super = 1
        return M_sub, M_super
    
    # Subsonic
    m = 0.5 # first guess
    for iter in range(1,101):
        eq = Aratio*m - ((1+(gam-1)/2*m**2)/((gam+1)/2))**((gam+1)/(gam-1)/2)  # Equation 9.44 in Cengal (without typo)
        deq = Aratio - (m*(((1/2*gam - 1/2)*m**2 + 1)/(1/2*gam + 1/2))**((1/2*(gam + 1))/(gam - 1) - 1)*(1/2*gam - 1/2)*(gam + 1))/((1/2*gam + 1/2)*(gam - 1))
        m = m - eq/deq
        if np.abs(eq)<1E-8:
            break
        if m < 0 or m > 1:
            # No good
            m = np.random(1)
    if iter >99:
        print('trouble converging for subsonic root with Aratio ' +str(Aratio))
        M_sub = np.NaN
    else:
        M_sub = m
    # Supersonic
    m = 2 # guess
    # Newton Rhapson
    for i in range(1,101):
        eq = Aratio*m - ((1+(gam-1)/2*m**2)/((gam+1)/2))**((gam+1)/(gam-1)/2)  # Equation 9.44 in Cengal (without typo)
        deq = Aratio - (m*(((1/2*gam - 1/2)*m**2 + 1)/(1/2*gam + 1/2))**((1/2*(gam + 1))/(gam - 1) - 1)*(1/2*gam - 1/2)*(gam + 1))/((1/2*gam + 1/2)*(gam - 1))
        m = m - eq/deq
        if np.abs(eq)<1E-8:
            break
        if m < 1 or m > 1E3:
            # No good
            m = np.random(1)*999+1
    if iter>99:
        print('trouble converging for supersonic root with Aratio ' +str(Aratio))
        M_super = np.NaN
    else:
        M_super = m
    return M_sub, M_super

def MforPratio(Pratio,gam):
    '''
    Calculates the mach number given a stagnation to static pressure ratio
    
    Args:
    Pratio (float): static pressure ratio
    gam (float): Specific Heat Ratio

    Returns:
    M (float): Mach Number
    '''
    

    M = np.sqrt((Pratio**((gam-1)/gam)-1)*2/(gam-1))
    return M

def MshockforPratio(Pratio,gam):
    '''
    Calculates the mach number given a stagnation to static pressure ratio
    
    Args:
    Pratio (float): static pressure ratio
    gam (float): Specific Heat Ratio
    Returns:
    M (float): Mach Number across shock
    ''' 
    M = np.sqrt((Pratio*(gam+1)+gam-1)/(2*gam))
    return M

def MshockforV2s_a1(V2s_a1,gam):
    '''
    This function calculates the mach number given 
    the ratio of the downstream velocity to the upstream speed of sound
    It is assumed that fluid upstream of the shock is stationary
    
    Args:
    V2s_a1 (float): Stream Ratio
    gam (float): Specific Heat Ratio
    '''
    M = (gam+1)/4*np.abs(V2s_a1)+np.sqrt(((gam+1)/4*V2s_a1)**2+1)
    return M

def obliqueDeltaforTheta(m1,theta,gamma):
    '''
    Finding flow turn angle 

    Args:
    m1 (float) = upstream mach number
    theta (float) = oblique shock angle relative to the incoming flow [RAD]
    gam = specific heat ratio

    Returns:
    delta (float): Turn angle
    '''
    tanDelta = (2*cot(theta)*(m1**2*np.sin(theta)**2-1))/(m1**2*(gamma+np.cos(2*theta))+2)
    delta = np.arctan(tanDelta)
    return delta


def obliqueThetasforDelta(m1,delta,gamma):
    '''
    Finds the strong and weak shong angles (theta) that satisfies the oblique shock relation
    that give you a specific turning angle

    Args:
     m1 (float) = upstream mach number
     delta (float): Turn angle
     gam = specific heat ratio

     Returns:
     theta_strong (float) = oblique strong shock angle relative to the incoming flow [RAD]
     theta_weak (float) = oblique weak shock angle relative to the incoming flow [RAD]

    '''
    tandelta = np.tan(delta)
    # Newton Iteration to solve equation for weak shock
    theta = np.sin(1/m1)
    for iter in range(1,101):
        Eq = (2*cot(theta)*(m1**2*np.sin(theta)**2-1))/(m1**2*(gamma*np.cos(2*theta))+2)
        Err = Eq -tandelta
        dErrdTheta = (4*m1**2*np.sin(2*theta)*cot(theta)*(m1**2*np.sin(theta)**2 - 1))/((np.cos(2*theta) + gamma)*m1**2 + 2)**2 - (2*(cot(theta)**2 + 1)*(m1**2*np.sin(theta)**2 - 1))/((np.cos(2*theta) + gamma)*m1**2 + 2) + (4*m1**2*np.cos(theta)*cot(theta)*np.sin(theta))/((np.cos(2*theta) + gamma)*m1**2 + 2)
        theta = theta - Err/dErrdTheta
        theta = np.abs(theta)
        if (np.abs(Err) < 1.0e-5):
            break
    if iter > 99 or theta < np.arcsin(1/m1) or theta > np.pi/2 :
        print('did not converge to weak shock solution, M = ' + str(m1) +' \delta = ' +str(delta*180/np.pi))
        theta_weak = np.NaN
    else:
        theta_weak = theta

    # Newton Iteration to solve equation for strong shock
    theta = np.pi/2
    for iter in range(1,101):
        Eq = (2*cot(theta)*(m1**2*np.sin(theta)**2-1))/(m1**2*(gamma+np.cos(2*theta))+2)
        Err = Eq -tandelta
        dErrdTheta = (4*m1**2*np.sin(2*theta)*cot(theta)*(m1**2*np.sin(theta)**2 - 1))/((np.cos(2*theta) + gamma)*m1**2 + 2)**2 - (2*(cot(theta)**2 + 1)*(m1**2*np.sin(theta)**2 - 1))/((np.cos(2*theta) + gamma)*m1**2 + 2) + (4*m1**2*np.cos(theta)*cot(theta)*np.sin(theta))/((np.cos(2*theta) + gamma)*m1**2 + 2)
        theta = theta - Err/dErrdTheta
        theta = abs(theta)
        if (abs(Err) < 1.0e-5):
            break
    if iter > 99 or theta < np.arcsin(1/m1) or theta > np.pi/2:
        print('did not converge to strong shock solution, M = ' +str(m1) +' \delta = ' +str(delta*180/np.pi))
        theta_strong = np.NaN
    else:
        theta_strong = theta

    if theta_weak > theta_strong:
        print('strong and weak are inverted?')
    return theta_weak,theta_strong

def prandtlMeyerMforNu(nu,gam):
    '''
    Use prandtlMeyer relation to get mach from nu
    Args:
    Nu (float) [Rad]: Prandlt-Meyer Angle
    gam (float): Specific Heat Ratio

    Returns:
    M (float): Mach Number
    '''
    const = np.sqrt((gam+1)/(gam-1))

    #This is the symbolic prandtly meyer function
    M = 1.5
    for iter in range(1,101):
        NuOfM = const*np.atan(np.sqrt(M**2-1)/const) -np.atan(np.sqrt(M**2-1))
        Residual = NuOfM - nu
        dNudM = M/((M**2 - 1)**(1/2)*(((M**2 - 1)*(gam - 1))/(gam + 1) + 1)) - 1/(M*(M**2 - 1)**(1/2))
        M = M - Residual/dNudM
        
        if (np.abs(Residual) < 1.0e-5):
            break
    return M

def prandtlMeyerNuforM(M,gam):
    '''
    Calculates Prandtl-Meyer function 

    Args: 
    M (float): Mach Number
    gam (float): Specific Heat Ratio

    Returns:
    Nu (float) [RAD]: Prandtl-Meyer Function
    '''
    const = np.sqrt((gam+1)/(gam-1))
    nu = const*np.arctan(np.sqrt(M**2-1)/const) -np.arctan(np.sqrt(M**2-1))

    return nu

def rayleighsMsforT0(T0_Tstar,gam):
    
    '''
    Calculates Subsonic and supersonic mach numbers given T0/T0* for Rayleigh flow

    Args:
    T0_Tstar (float): Temperature Ratio (Stagnation to throat)
    gam (float): Specific Heat Ratio
    
    Returns:
    M_sub (float): Subsonic Mach Num
    M_super (float): Supersonic Mach Num
    '''
    M_sub = (-(T0_Tstar*gam - gam + gam*(1 - T0_Tstar)**(1/2) + (1 - T0_Tstar)**(1/2) - 1)/(T0_Tstar*gam**2 - gam**2 + 1))**(1/2)
    M_super = ((gam - T0_Tstar*gam + gam*(1 - T0_Tstar)**(1/2) + (1 - T0_Tstar)**(1/2) + 1)/(T0_Tstar*gam**2 - gam**2 + 1))**(1/2)
    return M_sub, M_super

def rayleighsT0forM(M,gam):
    
    '''
    Calculates P/P*, T/T*, T0/T* and V/V* for Rayleigh flow given a Mach Number

    Args:
    M (float): Mach Num
    gam (float): Specific Heat Ratio
    
    Returns:
    T0_T0star (float): Temperature Ratio (Stagnation to throat)
    P_Pstar (float): Presure Ratio (Stagnation to throat)
    T_Tstar (float): Temperature Ratio (Dynamic to throat)
    V_Vstar (float): Velocity Ratio (Stagnation to throat)
    '''
    rayleigh_ratios = {
        'P_Pstar': (1+gam)/(1+gam*M**2),
        'T_Tstar': (1+gam)**2*M**2/(1+gam*M**2)**2
    }
    rayleigh_ratios['T0_T0star'] = rayleigh_ratios['T_Tstar']*(2+(gam-1)*M**2)/(gam+1)
    rayleigh_ratios['P0_P0star'] =rayleigh_ratios['P_Pstar']*((2+(gam-1)*M**2)/(gam+1))**(gam/(gam-1))
    rayleigh_ratios['V_Vstar'] = np.sqrt(rayleigh_ratios['T_Tstar'])*M
    return rayleigh_ratios

def shock(M1,gam):
    '''
    Calculates flow property values after a shock

    Args:
    M1 (float): Initial Mach Number
    gam (float): Specific Heat Ratio

    Returns:
    P2_P1 (float): Pressure Ratio
    T2_T1 (float): Temperature Ratio
    rho2_rho1 (float): Density Ratio
    P01_P02 (float): Stagnation Pressure Ratio
    Astar2_Astar1 (float): Imaginary Throat to Throat
    '''
    M2 = np.sqrt((1+(gam-1)/2*M1**2)/(gam*M1**2-(gam-1)/2))
    shock_ratios = {
        'P2_P1': (1+gam*M1**2)/(1+gam*M2**2),
        'T2_T1': (1+(gam-1)/2*M1**2)/(1+(gam-1)/2*M2**2),
    }
    shock_ratios['rho2_rho1'] = shock_ratios['P2_P1']/shock_ratios['T2_T1']
    shock_ratios['P02_P01'] = shock_ratios[('P2_P1')]/shock_ratios[('T2_T1')]**(gam/(gam-1))
    shock_ratios['V_Vstar'] = M2/M1*(shock_ratios['T2_T1)'])^((gam+1)/(2*(gam-1)));
    return shock_ratios
    





