%% Initial Concept Analysis Software - Nathan Rand
% 09/23/2023
close all;
clc;

%Figure formatting preferences
set(groot,'defaultAxesFontSize',14);
set(groot, 'DefaultLineLineWidth', 2);

%Define our inputs
stages = ["Air Breathing", "Rocket"];
Isp = [-1, 300];
alpha = [1/2, 1/4];
pi_e = [0.12, 0.12];
DDeF_ratio = [0, 0];
eta_o = [0.5, -1];
hpr = [120000e3, -1];
mp = 45.4;

%Run our analysis script
disp(vehicle_mass_analysis(stages, Isp, alpha, pi_e, DDeF_ratio, eta_o, hpr, mp));

function [mi] = vehicle_mass_analysis(stages, Isp, alpha, pi_e, DDeF_ratio, eta_o, hpr, mp)
    %Define basic givens for the conceptual problem
    ro = 6378e3; % (m) Radius of Earth
    ri = (6378+12.192)*1e3; % (m)
    rf = (6378+251.46)*1e3; % (m)
    go = 9.81; % (m/s^2)

    %Define memory space for our results to be stored for each stage
    lambda = zeros(length(stages));
    gamma = zeros(length(stages));

    %Iterate over all of our stages
    for i = 1:length(stages)
        %Compute results for each stage
        if upper(stages(i)) == "ROCKET"
            lambda(i) = sqrt(go*ro)/(go*Isp(i)*(1-DDeF_ratio(i)));
            gamma(i) = 1/(exp(-sqrt(alpha(i))*lambda(i))-pi_e(i));
        elseif upper(stages(i)) == "AIR BREATHING"
            lambda(i) = (go*ro*(1-(1/2)*(ro/rf)))/(eta_o(i)*hpr(i)*(1-DDeF_ratio(i)));
            gamma(i) = 1/(exp(-alpha(i)*lambda(i))-pi_e(i)); 
        else
            disp("Non-Supported Stage Identifier Passed. Please pass either ROCKET or AIR BREATHING keys.");
        end
    end

    %Determine cumulative results for the entire vehicle
    gamma_tot = prod(gamma);
    gamma_tot = gamma_tot(1);
    mi = gamma_tot*mp;
end