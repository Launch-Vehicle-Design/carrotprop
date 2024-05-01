%% Sizing and Heat Transfer for Optimization and Analysis of RDREs (SHOAR) - Nathan Rand
% 02/12/2024
close all;
clc;
clear;

%Figure formatting preferences
set(0,'DefaultLineLineWidth',1.75);
set(0,'DefaultLineMarkerSize',14);
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultFigureColor',[1,1,1]);
set(0,'DefaultTextFontSize',14);
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultTextFontName','Times-Roman');
set(0,'DefaultAxesFontName','Times-Roman');


%Define problem wide variables
H_release = 12192; % (m)
H_orbit = 251.46e3; % (m)
Weight = 2267.962*9.81; % (N)
g = 9.81; % m/s^2
Ru = 8314; % (J/kgK)
[Ta, a, Pa, rho] = atmosisa(H_release);

%Performance calculations based on the method outlined in the
%following paper: https://arc.aiaa.org/doi/epdf/10.2514/1.A34313
ceam_out = CEA('reac', 'name', 'C2H4', 'wt%',100.0, 'name','N2O4','wt%',100.0,'prob','det', 'equilibrium', 't,k',550.,'p,atm',15,'output','transport','end','screen');
tc = 120/1e6; % Detonation Cycle Time
aeat = 300;
t = 0:1e-6:tc;
gamma = ceam_out.output.burned.gamma;
P1 = ceam_out.output.unburned.pressure*1e5;
To = ceam_out.output.burned.temperature;
PR = ceam_out.output.p_ratio;
lambda = log(PR)/tc;
P = P1*PR*exp(-lambda*t);
T = To*(P./P1(1)).^((gamma-1)/gamma);
density = 1.263; % Density term defined by Huzel and Huang (density impulse)
c_star = sqrt(gamma*(Ru/ceam_out.output.burned.mw)*T)./(gamma*sqrt((2/(gamma+1))^((gamma+1)/(gamma-1))));
CF_spike = sqrt(((2*gamma^2)/(gamma-1))*((2/(gamma+1))^((gamma+1)/(gamma-1)))*(1-((1/NPR(aeat, gamma))^((gamma-1)/gamma)))+aeat*((1/NPR(aeat, gamma))-(Pa./P)));
Isp = CF_spike.*c_star./g;
Isp_vol = Isp.*density;
fprintf("Average Isp for RDRE Engine: %3.1f s\n", mean(Isp));
fprintf("Average Density Impulse (Ispd) is: %3.1f s\n\n", mean(Isp_vol));

%Produce related detonation plots
figure();
plot(t*1e6, P*145/(1e6));
xlabel("Time (microseconds)");
ylabel("Pressure (psia)");
title("Detonation Pressure v. Time");

figure();
plot(t*1e6, T*1.8);
xlabel("Time (microseconds)");
ylabel("Temperature (R)");
title("Detonation Temperature v. Time");

figure();
plot(t*1e6, c_star);
xlabel("Time (microseconds)");
ylabel("Characteristic Exhaust Velocity C* (m/s)");
title("Characteristic Exhaust Velocity v. Time");

figure();
plot(t*1e6, CF_spike);
xlabel("Time (microseconds)");
ylabel("Aerospike Thrust Coefficient (CF)");
title("Thrust Coefficient (Spike) v. Time");

figure();
plot(t*1e6, Isp);
xlabel("Time (microseconds)");
ylabel("Specific Impulse (s)");
title("Specific Impulse v. Time");

figure();
plot(t*1e6, Isp_vol);
xlabel("Time (microseconds)");
ylabel("Volumetric Specific Impulse (s)");
title("Volumetric Specific Impulse v. Time");


%% Define support functionality
%Returns the stagnation pressure ratio based on 
function npr = NPR(aeat, gamma)
    M = m_aas(aeat, gamma, 1);
    disp(M);
    npr = (1+((gamma-1)/2)*M.^2).^(gamma/(gamma-1));
end