%% Engine Sizing Determination and Optimization (ESDO) - Nathan Rand
% 09/30/2023
close all;
clc;
clear;

%Figure formatting preferences
set(groot,'defaultAxesFontSize',14);
set(groot, 'DefaultLineLineWidth', 2);

%Define problem wide variables
H_release = 12192; % (m)
H_orbit = 251.46e3; % (m)
Weight = 2267.962*9.81; % (N)

%Run CEA multiple times (currently the web CEARUN version) in order to
%determine rough lengths for our nozzles

%Liquid Engine (N2H4 and N2O4)
%Chamber Pressure: 5000 psia
Pc = 34.47e7; % (Pa)
Rt = 0.0254; % meters (1 inch)
Astar = pi*Rt^2; % Throat Area (m^2)
R_curve = 0.382*Rt; % meters
cea_data = readtable("data/CEA/Liquid/liquid_cea.txt");
cea_data.p = cea_data.p*1e5; % Convert to Pascals
cea_data.isp = cea_data.isp/9.81; % Convert to s
Thrust = cea_data.cf*Pc*Astar; % (N)
length_vals = (Rt*(sqrt(cea_data.aeat)-1)+R_curve*((1/cosd(15))-1))/tand(15);
%Determine area ratios that are optimal for specific Pa values
Isp_percent = 0.8:0.01:1;
epsilon = zeros([1,length(Isp_percent)]);
for i=1:length(Isp_percent)
    for j=2:length(cea_data.isp)
        disp(cea_data.isp(j))
        if(cea_data.isp(j) >= Isp_percent(i)*max(cea_data.isp))
            epsilon(i) = cea_data.aeat(j-1) + (Isp_percent(i)*max(cea_data.isp)-cea_data.isp(j-1))*((cea_data.aeat(j) - cea_data.aeat(j-1))/(cea_data.isp(j)-cea_data.isp(j-1)));
            break
        end
    end
end
disp(epsilon);
length_trunc = (Rt*(sqrt(epsilon)-1)+R_curve*((1/cosd(15))-1))/tand(15);

%Plot our results for length distributions
figure();
plot(cea_data.aeat, length_vals);
xlabel("Ae/At");
ylabel("Nozzle Length L_n (m)");
title("Nozzle Length v. Area Ratio");


%Plot our results for Isp vs area ratio
figure();
plot(cea_data.aeat, cea_data.isp);
xlabel("Ae/At");
ylabel("Specific Impulse Isp (s)");
title("Specific Impulse v. Area Ratio");

%Plot our results for Isp Percentage and Length
figure();
plot(Isp_percent*100, length_trunc);
xlabel("Specific Impulse Percentage");
ylabel("Nozzle Length");
title("Specific Impulse Percentage and Nozzle Length Truncation");

%Plot results for Thrust
figure();
plot(cea_data.aeat, Thrust/Weight);
xlabel("Ae/At");
ylabel("Thrust/Weight");
title("Thrust to Weight v. Area Ratio");
