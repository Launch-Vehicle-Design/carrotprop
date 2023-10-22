%% Engine Sizing Determination and Optimization (ESDO) - Nathan Rand
% 09/30/2023
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

%Define problem wide variables
H_release = 12192; % (m)
H_orbit = 251.46e3; % (m)
Weight = 2267.962*9.81; % (N)
g = 9.81; % m/s^2
Ru = 8314; % (J/kgK)
[Ta, a, Pa, rho] = atmosisa(H_release);

%Run CEA multiple times (currently the web CEARUN version) in order to
%determine rough lengths for our nozzles
mode = "RDRE";
cea_data = -1; % Default instantiation
if(mode == "File")
    cea_data = readtable("data/CEA/Liquid/rppero_data.txt");
    cea_data.isp = cea_data.isp/9.81; % Convert to s
    cea_data.ivac = cea_data.ivac/9.81; % Convert to s
else
    ceam_out = -1;
    if(mode == "CL20")
        ceam_out = CEA('reac','name','CL20','C', 6.28, 'H', 6.20, 'N', 12.55, 'O', 12, 'h,Kj/mol', 581.87,'wt%',100.0,'t(k)',298.15,'prob','rkt','p,psia',2900,'supar',5,10,15,20,25,30,35,40,45,50,55,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,'outp','massf','transport','mks','end');
        cea_data = array2table(squeeze([ceam_out.output.froz.aeat(3:length(ceam_out.output.froz.aeat)), ceam_out.output.froz.pressure(3:length(ceam_out.output.froz.pressure)), ...
                            ceam_out.output.froz.density(3:length(ceam_out.output.froz.density)), ceam_out.output.froz.mach(3:length(ceam_out.output.froz.mach)), ...
                            ceam_out.output.froz.sonvel(3:length(ceam_out.output.froz.sonvel)), ceam_out.output.froz.isp(3:length(ceam_out.output.froz.isp)), ...
                            ceam_out.output.froz.cf(3:length(ceam_out.output.froz.cf)), ceam_out.output.froz.isp_vac(3:length(ceam_out.output.froz.isp_vac))])');
        cea_data.Properties.VariableNames = ["aeat", "p", "rho", "mach", "son", "isp", "cf", "ivac"];
    elseif(mode == "RDRE")
        %Performance calculations based on the method outlined in the
        %following paper: https://arc.aiaa.org/doi/epdf/10.2514/1.A34313
        ceam_out = CEA('reac','name','H2','moles',8.000,'name','O2','moles',1.0000,'prob','det','t,k',500.,'p,atm',50.0,'output','transport','end','screen');
        tc = 120/1e6; % Detonation Cycle Time
        aeat = 5;
        t = 0:1e-6:tc;
        gamma = ceam_out.output.burned.gamma;
        P1 = ceam_out.output.unburned.pressure*1e5;
        To = ceam_out.output.burned.temperature;
        PR = ceam_out.output.p_ratio;
        lambda = log(PR)/tc;
        P = P1*PR*exp(-lambda*t);
        T = To*(P./P1(1)).^((gamma-1)/gamma);
        SG = 0.28; % From LVD Textbook for H2/O2
        c_star = sqrt(gamma*(Ru/ceam_out.output.burned.mw)*T)./(gamma*sqrt((2/(gamma+1))^((gamma+1)/(gamma-1))));
        CF_spike = sqrt(((2*gamma^2)/(gamma-1))*((2/(gamma+1))^((gamma+1)/(gamma-1)))*(1-((1/NPR(aeat, gamma))^((gamma-1)/gamma)))+aeat*((1/NPR(aeat, gamma))-(Pa./P)));
        Isp = CF_spike.*c_star./g;
        Isp_vol = Isp.*SG;

        %Produce related detonation plots
        figure();
        plot(t*1e6, P/1e6);
        xlabel("Time (microseconds)");
        ylabel("Pressure (MPa)");
        title("Detonation Pressure v. Time");

        figure();
        plot(t*1e6, T);
        xlabel("Time (microseconds)");
        ylabel("Temperature (K)");
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
    end
end
cea_data.p = cea_data.p*1e5; % Convert to Pascals

%Chamber Pressure: 10 MPa
Pc = 20e6; % (Pa)
Rt = 0.0254; % meters (3/4 inch)
Astar = pi*Rt^2; % Throat Area (m^2)
R_curve = 0.382*Rt; % meters
Thrust = cea_data.cf*Pc*Astar; % (N)
length_vals = (Rt*(sqrt(cea_data.aeat)-1)+R_curve*((1/cosd(15))-1))/tand(15);
%Determine area ratios that are optimal for specific Pa values
Isp_40 = (Astar*(cea_data.rho.*cea_data.aeat.*(cea_data.mach.*cea_data.son).^2 + cea_data.aeat.*(cea_data.p-Pa)))./(Astar.*cea_data.rho.*cea_data.aeat.*(cea_data.mach.*cea_data.son).*g);
Isp_eff = 0.95*(Isp_40 + (2/3)*(cea_data.ivac - Isp_40));
Isp_percent = 0.9:0.0001:1;
epsilon = zeros([1,length(Isp_percent)]);
for i=1:length(Isp_percent)
    for j=2:length(Isp_eff)
        if(Isp_eff(j) >= Isp_percent(i)*max(Isp_eff))
            epsilon(i) = cea_data.aeat(j-1) + (Isp_percent(i)*max(Isp_eff)-Isp_eff(j-1))*((cea_data.aeat(j) - cea_data.aeat(j-1))/(Isp_eff(j)-Isp_eff(j-1)));
            break
        end
    end
end
length_trunc = (Rt*(sqrt(epsilon)-1)+R_curve*((1/cosd(15))-1))/tand(15);

%Output primary results in console
fprintf("Nozzle Geometry Points -> Rt: %0.4f m | Re: %1.4f m\n", [Rt, sqrt(epsilon(find(Isp_percent==0.90))*Rt^2)]);
fprintf("Nozzle length for 90 Percent Efficiency: %0.3f m. Effective Isp at this value is: %0.3f s \n\n", [length_trunc(find(Isp_percent==0.90)), 0.90*max(Isp_eff)]);
fprintf("Nozzle Geometry Points -> Rt: %0.4f m | Re: %1.4f m\n", [Rt, sqrt(epsilon(find(Isp_percent==0.95))*Rt^2)]);
fprintf("Nozzle length for 95 Percent Efficiency: %0.3f m. Effective Isp at this value is: %0.3f s \n\n", [length_trunc(find(Isp_percent==0.95)), 0.95*max(Isp_eff)]);
fprintf("Nozzle Geometry Points -> Rt: %0.4f m | Re: %1.4f m\n", [Rt, sqrt(epsilon(find(Isp_percent==0.98))*Rt^2)]);
fprintf("Nozzle length for 98 Percent Efficiency: %0.3f m. Effective Isp at this value is: %0.3f s \n\n", [length_trunc(find(Isp_percent==0.98)), 0.98*max(Isp_eff)]);
fprintf("Nozzle Geometry Points -> Rt: %0.4f m | Re: %1.4f m\n", [Rt, sqrt(epsilon(find(Isp_percent==0.995))*Rt^2)]);
fprintf("Nozzle length for 99.5 Percent Efficiency: %0.3f m. Effective Isp at this value is: %0.3f s \n", [length_trunc(find(Isp_percent==0.995)), 0.995*max(Isp_eff)]);
fprintf("Nozzle Geometry Points -> Rt: %0.4f m | Re: %1.4f m\n", [Rt, sqrt(epsilon(find(Isp_percent==1))*Rt^2)]);
fprintf("Nozzle length for 100 Percent Efficiency: %0.3f m. Effective Isp at this value is: %0.3f s \n", [length_trunc(find(Isp_percent==1)), 1*max(Isp_eff)]);

%Plot our results for length distributions
figure();
plot(cea_data.aeat, length_vals);
xlabel("Ae/At");
ylabel("Nozzle Length $L_n$ (m)");
title("Nozzle Length v. Area Ratio");


%Plot our results for Isp vs area ratio
figure();
plot(cea_data.aeat, cea_data.isp);
hold on;
plot(cea_data.aeat, Isp_eff);
hold on;
plot(cea_data.aeat, cea_data.ivac);
xlabel("Ae/At");
ylabel("Specific Impulse Isp (s)");
title("Specific Impulse v. Area Ratio");
legend(["CEA SL", "Effective", "CEA Vacuum"]);

%Plot our results for Isp Percentage and Length
figure();
plot(Isp_percent*100, length_trunc);
xlabel("Effective Specific Impulse ($I_{sp, eff}$) Percentage");
ylabel("Nozzle Length $L_n$ (m)");
title("Specific Impulse Percentage and Nozzle Length Truncation");

%Plot results for Thrust
figure();
plot(cea_data.aeat, Thrust/Weight);
xlabel("Ae/At");
ylabel("Thrust/Weight");
title("Thrust to Weight v. Area Ratio");

%Write our results to an Excel Sheet for future reference
out_table = table(Isp_percent', length_trunc');
out_table.Properties.VariableNames = ["Isp_p", "Ln"];
writetable(out_table, 'out/solid/cl20_isp_length.xlsx');


%% Define support functionality
%Returns the stagnation pressure ratio based on 
function npr = NPR(aeat, gamma)
    M = m_aas(aeat, gamma, 1);
    npr = (1+((gamma-1)/2)*M^2)^(gamma/(gamma-1));
end