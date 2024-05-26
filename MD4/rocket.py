import numpy as np

def rocket_rhs(t, y, P0, Agrain, rhof, Tf, Rb, gamb, Astar, Ae, a, n, mr, g, v, Cd, Ar):
    x = y[0]
    v = y[1]
    mp = y[2]
    
    # Assuming StandardConditions function is defined
    T, Pb, rho = StandardConditions(x)
    
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
