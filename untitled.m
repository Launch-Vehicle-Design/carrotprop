%% Interfacial Fluid Mechanics - Nathan Rand
% 01/31/2023
close all;
clc;

%Figure formatting preferences
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultLineMarkerSize',14);
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultFigureColor',[1,1,1]);
set(0,'DefaultTextFontSize',14);
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultTextFontName','Times-Roman');
set(0,'DefaultAxesFontName','Times-Roman');

%% Question 2
clear;
fprintf("Question 2: \n-----------\n\n");

%Define givens and solve for each case

%Part A
Ro = 10e-6; % (m)
rho = 1000; % (kg/m^3)
gamma = 72.8e-3; % (N/m)
Ua = 0.2*sqrt(gamma/(rho*Ro)); % (m/s)
fprintf("(A) U: %0.4f m/s\n", Ua);

%Part B
Ro = 100e-6; % (m)
rho = 1000; % (kg/m^3)
gamma = 72.8e-3; % (N/m)
Ub = 0.2*sqrt(gamma/(rho*Ro)); % (m/s)
fprintf("(B) U: %0.4f m/s\n", Ub);

%Part C
Ro = 1e-3; % (m)
rho = 1000; % (kg/m^3)
gamma = 72.8e-3; % (N/m)
Uc = 0.2*sqrt(gamma/(rho*Ro)); % (m/s)
fprintf("(C) U: %0.4f m/s\n", Uc);

%Part D
Ro = 100e-6; % (m)
rho = 789; % (kg/m^3)
gamma = 22.4e-3; % (N/m)
Ud = 0.2*sqrt(gamma/(rho*Ro)); % (m/s)
fprintf("(D) U: %0.4f m/s\n", Ud);

%Part E
Ro = 100e-6; % (m)
rho = 13545.8; % (kg/m^3)
gamma = 486.5e-3; % (N/m)
Ue = 0.2*sqrt(gamma/(rho*Ro)); % (m/s)
fprintf("(E) U: %0.4f m/s\n", Ue);