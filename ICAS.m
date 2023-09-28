%% Initial Concept Analysis Software - Nathan Rand
% 09/23/2023
close all;
clc;

%Figure formatting preferences
set(groot,'defaultAxesFontSize',14);
set(groot, 'DefaultLineLineWidth', 2);

%Problem constants
v_vehicle = 7.43224; % m^3 (80 ft^2)

%Define our inputs for the Rocket -> Air Breathing -> Liquid Rocket Case
stages = ["Rocket", "Air Breathing", "Rocket"];
Isp = [300, -1, 425];
alpha = [1/9, 1/3, 1/9];
pi_e = [0.1, 0.1, 0.1];
DDeF_ratio = [0, 0, 0];
eta_o = [-1, 0.5, -1];
hpr = [-1, 48000e3, -1];
mp = 45.4;

%Run our analysis script
mi = vehicle_mass_analysis(stages, Isp, alpha, pi_e, DDeF_ratio, eta_o, hpr, mp);
fprintf("Initial Mass Required for Solid Rocket -> Air Breathing -> Liquid Rocket is: %4.0f kg \n", mi);

%Define our inputs for the Rocket -> Air Breathing -> Rocket Case
stages = ["Rocket", "Air Breathing", "Rocket"];
Isp = [300, -1, 300];
alpha = [1/9, 1/3, 1/9];
pi_e = [0.1, 0.1, 0.1];
DDeF_ratio = [0, 0, 0];
eta_o = [-1, 0.5, -1];
hpr = [-1, 48000e3, -1];
mp = 45.4;

%Run our analysis script
mi = vehicle_mass_analysis(stages, Isp, alpha, pi_e, DDeF_ratio, eta_o, hpr, mp);
fprintf("Initial Mass Required for Solid Rocket -> Air Breathing -> Solid Rocket is: %4.0f kg \n", mi);

%Define our inputs for the Rocket -> Rocket -> Rocket Case
stages = ["Rocket", "Rocket", "Rocket"];
Isp = [320, 320, 320];
alpha = [1/9, 1/9, 1/9];
pi_e = [0.1, 0.1, 0.1];
DDeF_ratio = [0, 0, 0];
eta_o = [-1, -1, -1];
hpr = [-1, -1, -1];
mp = 45.4;
rho = 2004; %CL20 Density (kg/m^3)

%Run our analysis script
mi = vehicle_mass_analysis(stages, Isp, alpha, pi_e, DDeF_ratio, eta_o, hpr, mp);
mf = mi - mi*pi_e(1) - mp;
v_required = mf/rho; % Required Volume
fprintf("Initial Mass Required for Solid Rocket -> Solid Rocket -> Solid Rocket is: %4.0f kg \n", mi);
fprintf("Required fuel mass (mf): %4.0f kg \n", mf);
fprintf("Required volume (V): %1.3f m^3 \n", v_required);
fprintf("Percent of Vehicle Volume: %2.2f percent \n\n", (v_required/v_vehicle)*100);

%Define our inputs for the Rocket -> Rocket -> Liquid
stages = ["Rocket", "Rocket", "Rocket"];
Isp = [300, 300, 425];
alpha = [1/9, 1/9, 1/9];
pi_e = [0.1, 0.1, 0.1];
DDeF_ratio = [0, 0, 0];
eta_o = [-1, -1, -1];
hpr = [-1, -1, -1];
mp = 45.4;

%Run our analysis script
mi = vehicle_mass_analysis(stages, Isp, alpha, pi_e, DDeF_ratio, eta_o, hpr, mp);
fprintf("Initial Mass Required for Solid Rocket -> Solid Rocket -> Liquid Rocket is: %4.0f kg \n", mi);

%Define our inputs for the Rocket -> Rocket Case
stages = ["Rocket", "Rocket"];
Isp = [300, 300];
alpha = [1/4, 1/4];
pi_e = [0.1, 0.1];
DDeF_ratio = [0, 0];
eta_o = [-1, -1];
hpr = [-1, -1];
mp = 45.4;

%Run our analysis script
mi = vehicle_mass_analysis(stages, Isp, alpha, pi_e, DDeF_ratio, eta_o, hpr, mp);
fprintf("Initial Mass Required for Solid Rocket -> Solid Rocket is: %4.0f kg \n", mi);

%Define our inputs for the Air Breathing -> Rocket Case
stages = ["Air Breathing", "Rocket"];
Isp = [-1, 300];
alpha = [1/2, 1/4];
pi_e = [0.1, 0.1];
DDeF_ratio = [0, 0];
eta_o = [0.5, -1];
hpr = [48000e3, -1];
mp = 45.4;

%Run our analysis script
mi = vehicle_mass_analysis(stages, Isp, alpha, pi_e, DDeF_ratio, eta_o, hpr, mp);
fprintf("Initial Mass Required for Airbreathing -> Solid Rocket is: %4.0f kg \n", mi);

%Define our inputs for the Air Breathing -> Liquid Rocket Case
stages = ["Air Breathing", "Rocket"];
Isp = [-1, 425];
alpha = [1/2, 1/4];
pi_e = [0.1, 0.1];
DDeF_ratio = [0, 0];
eta_o = [0.5, -1];
hpr = [120000e3, -1];
mp = 45.4;

%Run our analysis script
mi = vehicle_mass_analysis(stages, Isp, alpha, pi_e, DDeF_ratio, eta_o, hpr, mp);
fprintf("Initial Mass Required for Airbreathing -> Liquid Rocket is: %4.0f kg \n", mi);

function [mi] = vehicle_mass_analysis(stages, Isp, alpha, pi_e, DDeF_ratio, eta_o, hpr, mp)
    %Define basic givens for the conceptual problem
    ro = 6378e3; % (m) Radius of Earth
    ri = (6378+12.192)*1e3; % (m)
    rf = (6378+251.46)*1e3; % (m)
    go = 9.81; % (m/s^2)

    %Ensure that alpha values sum to 1
    alpha_sum = 0;
    for i = 1:length(stages)
        if upper(stages(i)) == "ROCKET"
            alpha_sum = alpha_sum + sqrt(alpha(i));
        elseif upper(stages(i)) == "AIR BREATHING"
            alpha_sum = alpha_sum + alpha(i);
        else
            error("Non-Supported Stage Identifier Passed. Please pass either ROCKET or AIR BREATHING keys.");
        end
    end
    if alpha_sum ~= 1
        error("Alpha values do not sum to 1 (considering for rocket they should be squared). Please check these inputs and try again.")
    end

    %Define memory space for our results to be stored for each stage
    lambda = zeros(length(stages));
    gamma = zeros(length(stages));

    %Iterate over all of our stages
    for i = 1:length(stages)
        %Compute results for each stage
        if upper(stages(i)) == "ROCKET"
            lambda(i) = (sqrt(go*ro))/(go*Isp(i)*(1-DDeF_ratio(i)));
            gamma(i) = 1/(exp(-sqrt(alpha(i))*lambda(i))-pi_e(i));
        elseif upper(stages(i)) == "AIR BREATHING"
            lambda(i) = (go*ro*(1-(1/2)*(ro/rf)))/(eta_o(i)*hpr(i)*(1-DDeF_ratio(i)));
            gamma(i) = 1/(exp(-alpha(i)*lambda(i))-pi_e(i)); 
        else
            error("Non-Supported Stage Identifier Passed. Please pass either ROCKET or AIR BREATHING keys.");
        end
    end

    %Determine cumulative results for the entire vehicle
    gamma_tot = prod(gamma);
    gamma_tot = gamma_tot(1);
    mi = gamma_tot*mp;
end