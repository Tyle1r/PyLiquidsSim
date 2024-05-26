import numpy as np
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
        dRdMp = Isp*g**2*(-1/mdot*xi+1/mdot +(Isp*xi-(Mr+Mp)/mdot)*1./(Mr+Mp))
        Mp = Mp -R/dRdMp

    return Mp

