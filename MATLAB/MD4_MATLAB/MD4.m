%% MD4
% Megan Albright 
% Sea Level
clc
clear all;
%Constants
mass_sat = 750;     %[kg]
alt = 500e3;        %[m]
acc_max = 12;       %[G's]
Qf = 4000e3;        %[J/kg]
rhof = 1800;        %[kg/m^3]
Tf = 4000;          %[K]
Rb = 0.400e3;       %[J/kgK]
gamb = 1.2;
n = 0.6;
a = 4/(100*1e6^n);  %[?]
nozz_th = 0.02;     %[m]
exp_an = 20;        %[deg]
rho_nozz = 800;     %[kg/m^3]
d = 4;              %[m] VARIABLE
height_max = 80;    %[m] VARIABLE
Cd = 0.4;
Pc = 16e6;      %[Pa] VARIABLE
Tp = 600;            %[seconds]
g = 9.81;            %[m/s^2]
[T,Pb,rhoatm] = StandardConditions(0);          %[Pa]
r_earth = 6.378e6;  %[m]
r_tot = r_earth + alt; %[m]
vel = sqrt(g*r_tot);      %[m/s]
P_atm = Pb;    %[Pa]
thick = 0.02;
Ar = pi/4*(d^2);
%% Also decide grain shape, area, throat and nozzle geometry, propellant needed
%Assume grain area constant 
%Chosen Constants
Pc_Patm = Pc/P_atm; 
mach = MforPratio(Pc_Patm,gamb);
ratiosA = isentropic(mach,gamb);
Ae_Astar = ratiosA.A_Astar;
% grain and throat area comes form Pc
ratiosforT = isentropic(1,gamb);
Tstar = Tf/ratiosforT.T0_T;             %[K]
cstar = sqrt(gamb*Rb*Tstar);            %[m/s]
Vstar = cstar*1;                        %[m/s]
Astar_Agrain = (((gamb+1)/2)^(gamb/(gamb-1))*Rb*Tstar*rhof*a*Pc^(n-1))/Vstar;
Astar = 0.01847;     %keep varying
Agrain = Astar / Astar_Agrain;
Aexit = Astar*Ae_Astar;
r_exit = sqrt(Aexit/pi);
h_nozz = r_exit / tand(exp_an);
vol_nozz = pi*r_exit*(r_exit+sqrt(h_nozz^2+r_exit^2))*thick;
mass_nozz = vol_nozz*rho_nozz;      %[kg]
mr = mass_nozz+mass_sat;            %[kg]

%% Use calc rocket to find Isp, mdot then use that to find Mp (eqn for centripetal)
[thrust,mdot,rdot] = calc_rocket(Pb,Agrain,rhof,Tf,Rb,gamb,Astar,Aexit,a,n);
Isp = thrust/(mdot*g);
E = 1/2*vel^2 + g*alt;
mp_est = calcM(Isp,g,mr,mdot,E);
mp = 1.635*mp_est;
m0 = mp+mr;
timeforfuel = mp/mdot;
weight = m0*g;
if weight > thrust
    disp('You aint goin anywhere')
else 
    disp ('We off the ground')
end

%% Calculations
odeopt = odeset('RelTol',1.0e-12,'AbsTol',1.0e-12,'Events',@myEventsFcn);
y0 = [0;0;mp];
[t,y] = ode45(@(t,y) rocket(t,y,Pc,Agrain,rhof,Tf,Rb,gamb,Astar,Aexit,a,n,mr,g,vel,Cd,Ar),[0 Tp],y0,odeopt);
plot(t,y(:,1))
legend("Ode45")
xlabel("Time")
ylabel("Altitude")
title("Height vs. Time")
%Validation - issues, mdot MATCH, ISP MATCH, rdot MATCH, 
% mp has sudden drop from 1.3e3 to 800 in one time step
x = g*Isp*(m0/mdot - (log(m0./(m0-mdot.*t))+1).*(m0/mdot-t)) - 1/2*g.*t.^2;
figure(2)
plot(t,y(:,3))
title("Propellant Mass vs. Time")
xlabel("Time")
ylabel("Propellant Mass")
figure(3)
plot(t,y(:,2))
hold on 
plot(t,vel*ones(size(t)))
title("Speed vs. Time")
xlabel("Time")
ylabel("Speed")
legend("Rocket Speed","Speed Needed to Stay in Orbit",'Location','southeast')

%% Finding Accel and G's
v = y(:,2);
for i = 2:length(t)-1
    accel(i,1) = (v(i+1)-v(i-1))/(t(i+1)-t(i-1));
end

for i = 1:length(t)-1
    tGs(i,1) = t(i,1);
end

figure(4)
plot(tGs,accel)
title("Acceleration vs. Time")
xlabel("Time")
ylabel("Acceleration")

figure(5)
Gs = accel./g;
plot(tGs,Gs)
title("G's vs. Time")
xlabel("Time")
ylabel("G's")
disp(max(Gs))
disp(max(v))
if max(v) > vel
    disp("We flyin")
else
    disp("we fallin back to the earf")
end

if max(Gs) > 12
    disp("satellite brokey")
else
    disp("satellite good :)")
end

%% Figure out height of rocket
vol_fuel = mp / rhof;
height = vol_fuel / (pi/4*d^2);



