# Adapted from Megan Albright's MATLAB Script
import utils.gasDynamics_utils
import MD4_utils.MD4_utils
from scipy.integrate import solve_ivp
import pandas as pd
import numpy as np

'''
This script takes a payload mass, the required altitude, and fuel properties
to simulate a rocket getting a payload to orbit
'''

mass_payload = 750     #[kg]
alt = 500e3        #[m]
max_accel = 12       #[G's]
Qf = 4000e3        #[J/kg]
rhof = 1800        #[kg/m**3]
Tf = 4000          #[K]
Rb = 0.400e3       #[J/kgK]
gamb = 1.2
n = 0.6
a = 4/(100*1e6**n)  #[m/s]
nozz_throatArea = 0.02     #[m**2]
exp_an = 20*(np.pi/180)        #[deg]
rho_nozz = 800     #[kg/m**3]
d = 4              #[m] VARIABLE
height_max = 80    #[m] VARIABLE
Cd = 0.4
Pc = 16e6      #[Pa] VARIABLE
Tp = 600            #[seconds]
g = 9.81            #[m/s**2]
T_atm,P_atm,rho_atm = utils.gasDynamics_utils.StandardConditions(0)          #[Pa]
r_earth = 6.378e6  #[m]
r_tot = r_earth + alt #[m]
vel = np.sqrt(g*r_tot)      #[m/s]
thick = 0.02
Ar = np.pi/4*(d**2)
## Also decide grain shape, area, throat and nozzle geometry, propellant needed
#Assume grain area constant 
#Chosen Constants
Pc_Patm = Pc/P_atm 
mach = utils.gasDynamics_utils.MforPratio(Pc_Patm,gamb)
ratiosA = utils.gasDynamics_utils.isentropic(mach,gamb)
Ae_Astar = ratiosA['A_Astar']
# grain and throat area comes form Pc
ratiosforT = utils.gasDynamics_utils.isentropic(1,gamb)
Tstar = Tf/ratiosforT['T0_T']             #[K]
cstar = np.sqrt(gamb*Rb*Tstar)            #[m/s]
Vstar = cstar*1                        #[m/s]
Astar_Agrain = (((gamb+1)/2)**(gamb/(gamb-1))*Rb*Tstar*rhof*a*Pc**(n-1))/Vstar
Astar = 0.01847     #keep varying
Agrain = Astar / Astar_Agrain
Aexit = Astar*Ae_Astar
r_exit = np.sqrt(Aexit/np.pi)
h_nozz = r_exit / np.tan(exp_an)
vol_nozz = np.pi*r_exit*(r_exit+np.sqrt(h_nozz**2+r_exit**2))*thick
mass_nozz = vol_nozz*rho_nozz      #[kg]
mr = mass_nozz+mass_payload            #[kg]

## Use calc rocket to find Isp, mdot then use that to find Mp (eqn for centripetal)
thrust,mdot,rdot = MD4_utils.MD4_utils.calc_rocket(P_atm,Agrain,rhof,Tf,Rb,gamb,Astar,Aexit,a,n)
Isp = thrust/(mdot*g)
E = 1/2*vel**2 + g*alt
mp_est = MD4_utils.MD4_utils.calcM(Isp,g,mr,mdot,E)
mp = 1.635*mp_est
m0 = mp+mr
timeforfuel = mp/mdot
weight = m0*g
if weight > thrust:
    print('You aint goin anywhere')
else: 
    print('We off the ground')

## Calculations
    '''
odeopt = {'rtol': 1.0e-12, 'atol': 1.0e-12, 'events': MD4_utils.MD4_utils.flightEvents}
'''
y0 = [0, 0, mp]

sol = solve_ivp(
    lambda t, y: MD4_utils.MD4_utils.rocket(t, y, Pc, Agrain, rhof, Tf, Rb, gamb, Astar, Aexit, a, n, mr, g, vel, Cd, Ar), 
    [0, Tp], y0, method='RK45', rtol=1.0e-12, atol=1.0e-12, events=lambda t, y: MD4_utils.MD4_utils.flightEvents(t, y, alt)
)
t = sol.t
y = sol.y
# Plot 



#Validation - issues, mdot MATCH, ISP MATCH, rdot MATCH, 
# mp has sudden drop from 1.3e3 to 800 in one time step
x = g*Isp*(m0/mdot - (np.log(m0/(m0-mdot*t))+1)*(m0/mdot-t)) - 1/2*g*t**2
'''
plot
title("Speed vs. Time")
xlabel("Time")
ylabel("Speed")
legend("Rocket Speed","Speed Needed to Stay in Orbit",'Location','southeast')
'''
## Finding Accel and G's
v = y[1,:]  # Assuming y is a 2D array
accel = np.zeros((len(t), 1))  # Initialize the acceleration array

for i in range(1, len(t) - 1):
    accel[i, 0] = (v[i + 1] - v[i - 1]) / (t[i + 1] - t[i - 1])
tGs = np.zeros((len(t), 1))  # Initialize the array
for i in range(1,len(t)-1):
    tGs[i,1] = t[i,1]
'''
figure(4)
plot(tGs,accel)
title("Acceleration vs. Time")
xlabel("Time")
ylabel("Acceleration")

figure(5)
Gs = accel./g
plot(tGs,Gs)
title("G's vs. Time")
xlabel("Time")
ylabel("G's")
print(max(Gs))
print(max(v))

'''
if max(v) > vel:
    print("We flyin")
else:
    print("we fallin back to the earf")

if max(tGs) > 12:
    print("satellite brokey")
else:
    print("satellite good :)")

## Figure out height of rocket
vol_fuel = mp / rhof
height = vol_fuel / (np.pi/4*d**2)