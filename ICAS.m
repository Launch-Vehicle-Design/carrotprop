%% Initial Concept Analysis Software - Nathan Rand
% 09/23/2023
close all;
clc;

%Figure formatting preferences
set(groot,'defaultAxesFontSize',14);
set(groot, 'DefaultLineLineWidth', 2);

%% Three Stage (Solid -> Solid -> Solid)
clear;
fprintf("Three Stage (Solid -> Solid -> Solid) Analysis: \n---------------------------------------------------------\n\n");

%Define basic givens for the conceptual problem
ro = 6378e3; % (m) Radius of Earth
ri = (6378+12.192)*1e3; % (m)
rf = (6378+251.46)*1e3; % (m)
go = 9.81; % (m/s^2)
DDeF_ratio = 0; % From (Heiser and Pratt) - rough approximation
Isp1 = 280; % All 3 stages (rough guess)
Isp2 = Isp1;
Isp3 = Isp1;
% Initial guess values for alphas and pi_e
alpha1 = 1/9;
alpha2 = 1/9;
alpha3 = 1/9;
pi_e = 0.12;
mp = 45.4; % (kg)

lambda1 = sqrt(go*ro)/(go*Isp1*(1-DDeF_ratio));
fprintf("Lambda 1 is: %1.3f \n", lambda1);
lambda2 = sqrt(go*ro)/(go*Isp2*(1-DDeF_ratio));
fprintf("Lambda 2 is: %1.3f \n", lambda2);
lambda3 = sqrt(go*ro)/(go*Isp3*(1-DDeF_ratio));
fprintf("Lambda 3 is: %1.3f \n", lambda3);

gamma1 = 1/(exp(-sqrt(alpha1)*lambda1)-pi_e);
fprintf("Gamma 1 is: %1.3f \n", gamma1);
gamma2 = 1/(exp(-sqrt(alpha2)*lambda2)-pi_e);
fprintf("Gamma 2 is: %1.3f \n", gamma2);
gamma3 = 1/(exp(-sqrt(alpha3)*lambda3)-pi_e);
fprintf("Gamma 3 is: %1.3f \n", gamma3);
gamma = gamma1*gamma2*gamma3;
fprintf("Overall Gamma: %3.1f \n", gamma);
mi1 = gamma*mp;
fprintf("Initial mass of the vehicle must be (mi1): %4.1f kg \n", mi1);

