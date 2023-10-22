%% Nozzle Length Processing Script (NLPS) - Nathan Rand
% 10/04/2023
close all;
clc;
clear;

%Figure formatting preferences
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultLineMarkerSize',14);
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultFigureColor',[1,1,1]);
set(0,'DefaultTextFontSize',14);
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultTextFontName','Times-Roman');
set(0,'DefaultAxesFontName','Times-Roman');

%Read in our data sets
hydT = readtable('out/liquid/hydnit_isp_length.xlsx');
hydT.Properties.VariableNames = ["Isp_percent", "Ln"];
rpT = readtable('out/liquid/rppero_isp_length.xlsx');
rpT.Properties.VariableNames = ["Isp_percent", "Ln"];
clT = readtable('out/solid/cl20_isp_length.xlsx');
clT.Properties.VariableNames = ["Isp_percent", "Ln"];

%Produce plot compare isp length percentages
figure();
plot(hydT.Isp_percent, hydT.Ln);
hold on;
plot(rpT.Isp_percent, rpT.Ln);
hold on;
plot(clT.Isp_percent, clT.Ln);
xlabel("Effective Specific Impulse ($I_{sp, eff}$) Percentage");
ylabel("Nozzle Length $L_n$ (m)");
title("Specific Impulse Percentage and Nozzle Length Truncation");
legend(["N2H4/N2O4 (RDRE)", "RP1/H2O2 (RDRE)", "CL20 (Solid)"]);