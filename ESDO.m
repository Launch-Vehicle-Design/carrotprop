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
mode = "CL20";
binder = "GAP";
cea_data = -1; % Default instantiation
if(mode == "File")
    cea_data = readtable("data/CEA/Liquid/rppero_data.txt");
    cea_data.isp = cea_data.isp/9.81; % Convert to s
    cea_data.ivac = cea_data.ivac/9.81; % Convert to s
else
    ceam_out = -1;
    if(mode == "CL20")
        if binder == "HTPB"
            ceam_out = CEA('reac','name','CL20','C', 6.28, 'H', 6.20, 'N', 12.55, 'O', 12, 'h,Kj/mol', 581.87,'wt%', 88.0, 't(k)',298.15, 'name','HTPB','C', 213.8, 'H', 323, 'O', 4.6, 'N', 2.3, 'h,Kj/mol', 342.0, 'wt%', 12.0, 't(k)', 298.15, 'prob','rkt','p,psia',2900,'supar',5,10,15,20,25,30,35,40,45,50,55,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,300,'outp','massf','transport','mks','end');
        elseif binder == "GAP"
            ceam_out = CEA('reac','name','CL20','C', 6.28, 'H', 6.20, 'N', 12.55, 'O', 12, 'h,Kj/mol', 581.87,'wt%', 88.0, 't(k)',298.15, 'name','GAP','C', 388, 'H', 330, 'O', 4, 'N', 192, 'h,Kj/mol', 6060.0, 'wt%', 12.0, 't(k)', 298.15, 'prob','rkt','p,psia',2900,'supar',5,10,15,20,25,30,35,40,45,50,55,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,300,'outp','massf','transport','mks','end');
        end
        cea_data = array2table(squeeze([ceam_out.output.froz.aeat(3:length(ceam_out.output.froz.aeat)), ceam_out.output.froz.pressure(3:length(ceam_out.output.froz.pressure)), ...
                            ceam_out.output.froz.density(3:length(ceam_out.output.froz.density)), ceam_out.output.froz.mach(3:length(ceam_out.output.froz.mach)), ...
                            ceam_out.output.froz.sonvel(3:length(ceam_out.output.froz.sonvel)), ceam_out.output.froz.isp(3:length(ceam_out.output.froz.isp)), ...
                            ceam_out.output.froz.cf(3:length(ceam_out.output.froz.cf)), ceam_out.output.froz.isp_vac(3:length(ceam_out.output.froz.isp_vac))])');
        cea_data.Properties.VariableNames = ["aeat", "p", "rho", "mach", "son", "isp", "cf", "ivac"];
        cea_data.p = cea_data.p*1e5; % Convert to Pascals

        %Chamber Pressure: 10 MPa
        Pc = 20e6; % (Pa)
        Rt = 0.0185; % meters (approx. 0.72 inches)
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
        fprintf("Nozzle length for 99.5 Percent Efficiency: %0.3f m. Effective Isp at this value is: %0.3f s \n\n", [length_trunc(find(Isp_percent==0.995)), 0.995*max(Isp_eff)]);
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
    elseif(mode == "RDRE")
        %Performance calculations based on the method outlined in the
        %following paper: https://arc.aiaa.org/doi/epdf/10.2514/1.A34313
        ceam_out = CEA('reac', 'name', 'C2H4', 'wt%',100.0,'name','N2O4','wt%',100.0,'prob','det', 'equilibrium', 't,k',550.,'p,atm',50.0,'output','transport','end','screen');
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
    elseif mode == "Air-Breathing"
        %Scramjet Analysis
        %Assume we are flying at some Mach Number (M)
        Mo = 6;
        gamma = 1.4;
        R = 287;  % (J/kgK)
        Po = Pa;
        To = Ta;

        %Inlet calculations
        fprintf("Inlet conditions (Compression):\n--------------------------\n\n");

        %Define standards for our inlet behavior
        theta = 10; % (deg) (four 10 degree turns in inlet)
        
        %Across Shock 1
        fprintf("Across Shock 1: \n\n");
        %First determine our betaA value
        syms betaA
        func = 2*(1/tand(betaA))*((Mo^2*(sind(betaA)^2)-1)/(Mo^2*(gamma+cosd(2*betaA))+2)) - tand(theta);
        betaA = vpasolve(func, betaA, theta); % (deg)
        fprintf("Beta angle (B_A): %2.2f deg\n", betaA);
        %Now we compute values across the shock
        Mno = Mo*sind(betaA);
        fprintf("Mno: %1.3f\n", Mno);
        PA = Po*(1+(2*gamma/(gamma+1))*(Mno^2-1));
        fprintf("PA: %1.3f kPa \n", PA/1e3);
        Tto = To*(1+((gamma-1)/2)*Mo^2);
        fprintf("Tto: %4.0f K\n", Tto);
        Pto = Po*(1+((gamma-1)/2)*Mo^2)^(gamma/(gamma-1));
        fprintf("Pto: %1.3f MPa \n", Pto/1e6);
        PtA = Pto*(((1+(2*gamma/(gamma+1))*(Mno^2-1))^(-1/(gamma-1)))*(((2+(gamma-1)*Mno^2)/((gamma+1)*Mno^2))^(-gamma/(gamma-1))));
        fprintf("PtA: %1.3f MPa \n", PtA/1e6);
        TtA = Tto; % (K)
        MnA = sqrt((Mno^2 + (2/(gamma-1)))/((2*gamma/(gamma-1))*Mno^2-1));
        fprintf("MnA: %1.3f \n", MnA);
        MA = MnA/sind(betaA-theta);
        fprintf("MA: %1.3f \n", MA);
        TA = TtA*(1+((gamma-1)/2)*MA^2)^-1;
        fprintf("TA: %3.1f K \n\n", TA);
        
        %Across Shock 2
        fprintf("Across Shock 2: \n\n");
        %First determine our betaB value
        syms betaB
        func = 2*(1/tand(betaB))*((MA^2*(sind(betaB)^2)-1)/(MA^2*(gamma+cosd(2*betaB))+2)) - tand(theta);
        betaB = vpasolve(func, betaB, theta); % (deg)
        fprintf("Beta angle (B_B): %2.2f deg\n", betaB);
        %Now we compute values across the shock
        MnA = MA*sind(betaB);
        fprintf("MnA: %1.3f\n", MnA);
        PB = PA*(1+(2*gamma/(gamma+1))*(MnA^2-1));
        fprintf("PB: %1.3f kPa \n", PB/1e3);
        PtB = PtA*(((1+(2*gamma/(gamma+1))*(MnA^2-1))^(-1/(gamma-1)))*(((2+(gamma-1)*MnA^2)/((gamma+1)*MnA^2))^(-gamma/(gamma-1))));
        fprintf("PtB: %1.3f MPa \n", PtB/1e6);
        TtB = TtA; % (K)
        MnB = sqrt((MnA^2 + (2/(gamma-1)))/((2*gamma/(gamma-1))*MnA^2-1));
        fprintf("MnB: %1.3f \n", MnB);
        MB = MnB/sind(betaB-theta);
        fprintf("MB: %1.3f \n", MB);
        TB = TtB*(1+((gamma-1)/2)*MB^2)^-1;
        fprintf("TB: %3.1f K \n\n", TB);
        
        %Across Shock 3
        fprintf("Across Shock 3: \n\n");
        %First determine our betaC value
        syms betaC
        func = 2*(1/tand(betaC))*((MB^2*(sind(betaC)^2)-1)/(MB^2*(gamma+cosd(2*betaC))+2)) - tand(theta);
        betaC = vpasolve(func, betaC, theta); % (deg)
        fprintf("Beta angle (B_C): %2.2f deg\n", betaC);
        %Now we compute values across the shock
        MnB = MB*sind(betaC);
        fprintf("MnB: %1.3f\n", MnB);
        PC = PB*(1+(2*gamma/(gamma+1))*(MnB^2-1));
        fprintf("PC: %1.3f kPa \n", PC/1e3);
        PtC = PtB*(((1+(2*gamma/(gamma+1))*(MnB^2-1))^(-1/(gamma-1)))*(((2+(gamma-1)*MnB^2)/((gamma+1)*MnB^2))^(-gamma/(gamma-1))));
        fprintf("PtC: %1.3f MPa \n", PtC/1e6);
        TtC = TtB; % (K)
        MnC = sqrt((MnB^2 + (2/(gamma-1)))/((2*gamma/(gamma-1))*MnB^2-1));
        fprintf("MnC: %1.3f \n", MnC);
        MC = MnC/sind(betaC-theta);
        fprintf("MC: %1.3f \n", MC);
        TC = TtC*(1+((gamma-1)/2)*MC^2)^-1;
        fprintf("TC: %3.1f K \n\n", TC);
        
        %Across Shock 4
        fprintf("Across Shock 4: \n\n");
        %First determine our betaC value
        syms betaD
        func = 2*(1/tand(betaD))*((MC^2*(sind(betaD)^2)-1)/(MC^2*(gamma+cosd(2*betaD))+2)) - tand(theta);
        betaD = vpasolve(func, betaD, theta); % (deg)
        fprintf("Beta angle (B_D): %2.2f deg\n", betaD);
        %Now we compute values across the shock
        MnC = MC*sind(betaD);
        fprintf("MnC: %1.3f\n", MnC);
        PD = PC*(1+(2*gamma/(gamma+1))*(MnC^2-1));
        fprintf("PD: %1.3f kPa \n", PD/1e3);
        PtD = PtC*(((1+(2*gamma/(gamma+1))*(MnC^2-1))^(-1/(gamma-1)))*(((2+(gamma-1)*MnC^2)/((gamma+1)*MnC^2))^(-gamma/(gamma-1))));
        fprintf("PtD: %1.3f MPa \n", PtD/1e6);
        TtD = TtC; % (K)
        MnD = sqrt((MnC^2 + (2/(gamma-1)))/((2*gamma/(gamma-1))*MnC^2-1));
        fprintf("MnD: %1.3f \n", MnD);
        MD = MnD/sind(betaD-theta);
        fprintf("MD: %1.3f \n", MD);
        TD = TtD*(1+((gamma-1)/2)*MD^2)^-1;
        fprintf("TD: %3.1f K \n\n", TD);
        
        %Final results of the inlet compression
        psi = TD/To;
        fprintf("psi: %1.3f \n", psi);
        Press_ratio = (PD/Po)*((1/psi)^(gamma/(gamma-1)));
        fprintf("Pressure ratio (pi_c) is: %0.4f \n", Press_ratio);
        eta_c = (psi - (1/Press_ratio)^((gamma-1)/gamma))/(psi-1);
        fprintf("Compression Efficiency (eta_c): %0.4f \n", eta_c);
        P3 = double(PD);
        fprintf("P3: %1.3f kPa \n", P3/1e3);
        T3 = double(TD);
        fprintf("T3: %3.1f K \n", T3);
        M3 = MD;
        fprintf("M3: %1.3f \n", M3);
        u3 = M3*sqrt(gamma*R*T3);
        fprintf("u3: %3.1f m/s\n\n", u3);

        %Constant Pressure Heat Addition (assume u3 = u4)
        fprintf("Constant Pressure Combustion:\n--------------------------\n\n");
        ceam_out = CEA('reac','oxid','Air','wtfrac',1,'t(k)',T3,'fuel','RP-1','wtfrac', 1,'t(k)',298.15,'prob','case','Air-Breathing','hp','p(bar)', P3/1e5, 'o/f',14.72,'output','trace',1e-15,'mks','end','screen');
        Pt4 = PtD;
        u4 = u3;
        f = 1/ceam_out.output.oxfl;
        T4 = ceam_out.output.temperature;
        gamma = ceam_out.output.gamma;
        Cp = ceam_out.output.cp*1000;
        R = 8314/ceam_out.output.mw;

        %Flow Expansion (Nozzle)
        fprintf("\nFlow Expansion (Nozzle):\n--------------------------\n\n");
        M4 = u4/sqrt(gamma*R*T4);
        fprintf("M4: %1.3f \n", M4);
        Tt4 = T4*(1+((gamma-1)/2)*M4^2);
        fprintf("Tt4: %3.1f K \n", Tt4);
        T10 = Tt4*((Po/Pt4)^((gamma-1)/gamma));
        fprintf("T10: %3.1f K \n", T10);
        ue = sqrt(2*Cp*gamma*(Tt4-T10));
        fprintf("ue: %3.1f m/s \n", ue);
        ua = Mo*sqrt(1.4*287*To);
        fprintf("ua: %3.1f m/s \n", ua);
        specific_thrust = (1+f)*ue-ua;
        fprintf("Specific Thrust: %1.3f kN*s/kg \n", specific_thrust/1e3);
        Isp = specific_thrust/(f*g);
        fprintf("Isp: %3.1f s \n\n", Isp);

    end
end

%% Define support functionality
%Returns the stagnation pressure ratio based on 
function npr = NPR(aeat, gamma)
    M = m_aas(aeat, gamma, 1);
    disp(M);
    npr = (1+((gamma-1)/2)*M.^2).^(gamma/(gamma-1));
end