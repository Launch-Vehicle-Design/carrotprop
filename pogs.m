%% Propulsion Optimization and Geometry Software (POGS) - Nathan Rand
% 02/15/2024
close all;
clc;
clear;

%Figure formatting preferences
set(0,'DefaultLineLineWidth',2);
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
Weight = 1801.669*9.81; % (N)
g = 9.81; % m/s^2
Ru = 8314; % (J/kgK)
[Ta, a, Pa, rho] = atmosisa(H_release);

%Define run conditions
fuel = "CL20";
ox = "N2O4";
mode = "RDRE";
binder = "GAP";
ceam_out = -1;

%Run CEA in different configurations depending on the mode
if(mode == "Solid")
    mode = "analyze";
    if mode == "optimize"
        if fuel == "COTS"
            ceam_out = CEA('')
        end
        elseif fuel == "CL20"
            if binder == "HTPB"
                ceam_out = CEA('reac','name','CL20','C', 6.28, 'H', 6.20, 'N', 12.55, 'O', 12, 'h,Kj/mol', 581.87,'wt%', 88.0, 't(k)',298.15, 'name','HTPB','C', 213.8, 'H', 323, 'O', 4.6, 'N', 2.3, 'h,Kj/mol', 342.0, 'wt%', 12.0, 't(k)', 298.15, 'prob','rkt','p,psia',2900,'supar',5,10,15,20,25,30,35,40,45,50,55,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,300,'outp','massf','transport','mks','end');
            elseif binder == "GAP"
                ceam_out = CEA('reac','name','CL20','C', 6.28, 'H', 6.20, 'N', 12.55, 'O', 12, 'h,Kj/mol', 581.87,'wt%', 88.0, 't(k)',298.15, 'name','GAP','C', 388, 'H', 330, 'O', 4, 'N', 192, 'h,Kj/mol', 6060.0, 'wt%', 12.0, 't(k)', 298.15, 'prob','rkt','p,psia',2900,'supar',5,10,15,20,25,30,35,40,45,50,55,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,300,'outp','massf','transport','mks','end');
                %ceam_out = CEA('reac','name','CL20','C', 6.28, 'H', 6.20, 'N', 12.55, 'O', 12, 'h,Kj/mol', 581.87,'wt%', 88.0, 't(k)',298.15, 'name','GAP','C', 388, 'H', 330, 'O', 4, 'N', 192, 'h,Kj/mol', 6060.0, 'wt%', 12.0, 't(k)', 298.15, 'prob','rkt','p,psia',2900,'subar', 8, 6, 4, 2, 1.5, 1.25, 'supar',1.25, 1.5, 2, 2.5, 3.5, 5,10,12, 14.06,'outp','massf','transport','mks','end');
            end
        end
        cea_data = array2table(squeeze([ceam_out.output.froz.aeat(3:length(ceam_out.output.froz.aeat)), ceam_out.output.froz.pressure(3:length(ceam_out.output.froz.pressure)), ...
                            ceam_out.output.froz.density(3:length(ceam_out.output.froz.density)), ceam_out.output.froz.temperature(3:length(ceam_out.output.froz.temperature)), ceam_out.output.froz.mach(3:length(ceam_out.output.froz.mach)), ...
                            ceam_out.output.froz.sonvel(3:length(ceam_out.output.froz.sonvel)), ceam_out.output.froz.isp(3:length(ceam_out.output.froz.isp)), ...
                            ceam_out.output.froz.cf(3:length(ceam_out.output.froz.cf)), ceam_out.output.froz.isp_vac(3:length(ceam_out.output.froz.isp_vac))])');
        cea_data.Properties.VariableNames = ["aeat", "p", "rho", "T", "mach", "son", "isp", "cf", "ivac"];
        cea_data.p = cea_data.p*1e5; % Convert to Pascals
    
        %Chamber Pressure: 7.5 MPa
        Pc = 7.5e6; % (Pa)
        Rt = 0.07; % meters (approx. 0.72 inches)
        rho = 1948.8; % (kg/m^3)
        rb = 23.5/1e3; % (m/s)
        Astar = pi*Rt^2; % Throat Area (m^2)
        R_curve = 0.382*Rt; % meters
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
    
        %Calculate mass flow and thrust
        mdot = rb*rho*pi*(0.5588^2)/4;
        Thrust = mdot*Isp_eff*9.81;
        boost_mdot = mdot*1.2;
        boost_thrust = boost_mdot*Isp_eff*9.81;
        
        %Output primary results in console
        fprintf("Nozzle Geometry Points -> Rt: %0.4f m | Re: %1.4f | Ae/At: %2.2f\n", [Rt, sqrt(epsilon(find(Isp_percent==0.90))*Rt^2), epsilon(find(Isp_percent==0.90))]);
        fprintf("Nozzle length for 90 Percent Efficiency: %0.3f m. Effective Isp at this value is: %0.3f s \n\n", [length_trunc(find(Isp_percent==0.90)), 0.90*max(Isp_eff)]);
        fprintf("Nozzle Geometry Points -> Rt: %0.4f m | Re: %1.4f | Ae/At: %2.2f\n", [Rt, sqrt(epsilon(find(Isp_percent==0.95))*Rt^2), epsilon(find(Isp_percent==0.95))]);
        fprintf("Nozzle length for 95 Percent Efficiency: %0.3f m. Effective Isp at this value is: %0.3f s \n\n", [length_trunc(find(Isp_percent==0.95)), 0.95*max(Isp_eff)]);
        fprintf("Nozzle Geometry Points -> Rt: %0.4f m | Re: %1.4f | Ae/At: %2.2f\n", [Rt, sqrt(epsilon(find(Isp_percent==0.98))*Rt^2), epsilon(find(Isp_percent==0.98))]);
        fprintf("Nozzle length for 98 Percent Efficiency: %0.3f m. Effective Isp at this value is: %0.3f s \n\n", [length_trunc(find(Isp_percent==0.98)), 0.98*max(Isp_eff)]);
        fprintf("Nozzle Geometry Points -> Rt: %0.4f m | Re: %1.4f | Ae/At: %2.2f\n", [Rt, sqrt(epsilon(find(Isp_percent==0.995))*Rt^2), epsilon(find(Isp_percent==0.995))]);
        fprintf("Nozzle length for 99.5 Percent Efficiency: %0.3f m. Effective Isp at this value is: %0.3f s \n\n", [length_trunc(find(Isp_percent==0.995)), 0.995*max(Isp_eff)]);
        fprintf("Nozzle Geometry Points -> Rt: %0.4f m | Re: %1.4f | Ae/At: %2.2f\n", [Rt, sqrt(epsilon(find(Isp_percent==1))*Rt^2), epsilon(find(Isp_percent==1))]);
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
        hold on;
        plot(cea_data.aeat, boost_thrust/Weight);
        xlabel("Ae/At");
        ylabel("Thrust/Weight");
        title("Thrust to Weight v. Area Ratio");
        legend("Nominal Thrust", "Boost Thrust");
    
        %Plot physical properties (P, T, rho) v. axial position
        figure();
        plot(cea_data.aeat, cea_data.p/Pc);
        hold on;
        plot(cea_data.aeat, cea_data.T/cea_data.T(1));
        xlabel("Expansion Ratio (Ae/At)");
        ylabel("Stagnation Ratios");
        title("Stagnation Ratios v. Ae/At");
        legend("P/Pc", "T/Tc");
    elseif mode == "analyze"
        %Define basic geometry dimensions based on previous analysis
        Rc = 11.8; % (in)
        Rt = 2.756; % (in)
        Re = sqrt(15.25*(Rt^2)); % (in)
        Ln = 0.7263*12; % (in)
        n = 50; % Number of point to evaluate CEA at

        %Define our control points
        cp = [-(Rc-Rt)/tand(40), Rc; ...
                          0, Rt; ...
                          0.382*Rt*sind(18), 0.382*Rt*(1-cosd(18))+Rt;...
                          Ln, Re];
        xq = [linspace(cp(1,1), cp(2,1), n/2), linspace(cp(2,1), cp(4,1), n/2)];
        xq = unique(xq);
        pp = csape(cp(:,1)', [0, cp(:,2)', tand(8)], [1 1]);
        s = fnval(pp, xq);

        %Plot our generated TCA geometry [Thrust Chamber Assembly]
        figure();
        plot(xq, s, 'k-');
        hold on;
        plot(cp(:,1), cp(:,2), 'ro');
        xlabel("X (in)");
        ylabel("Y (in)");
        title("Solid Engine Contour");
        legend(["Cubic Spline Geometry", "Control Points"]);
        grid on;

        %Compute area ratio arrays for downstream and upstream of the
        %throat
        epsilon_down = [];
        epsilon_up = [];
        for i = 1:length(s)
            if xq(i) < 0
                epsilon_down = [epsilon_down, (s(i)^2)/(Rt^2)];
            elseif xq(i) >= 0
                epsilon_up = [epsilon_up, (s(i)^2)/(Rt^2)];
            end
        end

        %Remove high contraction ratio epsilons to avoid CEA failing
        epsilon_subsonic_limit = epsilon_down >= 8;
        epsilon_down(epsilon_subsonic_limit) = [];

        %Fit a cubic spline to the control points
        %Run CEA Analysis on the given geometry for the subsonic points
        for i = 1:length(epsilon_down)
            ceam_out = CEA('reac','name','CL20','C', 6.28, 'H', 6.20, 'N', 12.55, 'O', 12, 'h,Kj/mol', 581.87,'wt%', 21.80, 't(k)',298.15, 'name','GAP','C', 388, 'H', 330, 'O', 4, 'N', 192, 'h,Kj/mol', 6060.0, 'wt%', 12.0, 't(k)', 298.15, 'name', 'NH4CLO4(I)', 't(k)', 298.15, 'wt(%)', 43.20, 'name', 'AL', 't(k)', 298.15, 'wt(%)', 18.00, 'name', 'C2CL6', 't(k)', 298.15, 'wt(%)', 5.00,  'prob','rkt','p,psia',1088,'subar', epsilon_down(i),'outp','massf','transport','mks','end');
            if i == 1
                cea_data = array2table([ceam_out.output.froz.aeat(3:length(ceam_out.output.froz.aeat)), ceam_out.output.froz.pressure(3:length(ceam_out.output.froz.pressure)), ...
                            ceam_out.output.froz.density(3:length(ceam_out.output.froz.density)), ceam_out.output.froz.temperature(3:length(ceam_out.output.froz.temperature)), ceam_out.output.froz.mach(3:length(ceam_out.output.froz.mach)), ...
                            ceam_out.output.froz.sonvel(3:length(ceam_out.output.froz.sonvel)), ceam_out.output.froz.isp(3:length(ceam_out.output.froz.isp)), ...
                            ceam_out.output.froz.cf(3:length(ceam_out.output.froz.cf)), ceam_out.output.froz.isp_vac(3:length(ceam_out.output.froz.isp_vac))]);
                cea_data.Properties.VariableNames = ["aeat", "p", "rho", "T", "mach", "son", "isp", "cf", "ivac"];
            else
                cea_data_iter = array2table([ceam_out.output.froz.aeat(3:length(ceam_out.output.froz.aeat)), ceam_out.output.froz.pressure(3:length(ceam_out.output.froz.pressure)), ...
                            ceam_out.output.froz.density(3:length(ceam_out.output.froz.density)), ceam_out.output.froz.temperature(3:length(ceam_out.output.froz.temperature)), ceam_out.output.froz.mach(3:length(ceam_out.output.froz.mach)), ...
                            ceam_out.output.froz.sonvel(3:length(ceam_out.output.froz.sonvel)), ceam_out.output.froz.isp(3:length(ceam_out.output.froz.isp)), ...
                            ceam_out.output.froz.cf(3:length(ceam_out.output.froz.cf)), ceam_out.output.froz.isp_vac(3:length(ceam_out.output.froz.isp_vac))]);
                cea_data_iter.Properties.VariableNames = ["aeat", "p", "rho", "T", "mach", "son", "isp", "cf", "ivac"];
                cea_data = [cea_data; cea_data_iter];
            end
            
        end
        %Now do the same for the supersonic points
        for i = 1:length(epsilon_up)
            ceam_out = CEA('reac','name','CL20','C', 6.28, 'H', 6.20, 'N', 12.55, 'O', 12, 'h,Kj/mol', 581.87,'wt%', 21.80, 't(k)',298.15, 'name','GAP','C', 388, 'H', 330, 'O', 4, 'N', 192, 'h,Kj/mol', 6060.0, 'wt%', 12.0, 't(k)', 298.15, 'name', 'NH4CLO4(I)', 't(k)', 298.15, 'wt(%)', 43.20, 'name', 'AL', 't(k)', 298.15, 'wt(%)', 18.00,'C2CL6', 't(k)', 298.15, 'wt(%)', 5.0, 'prob','rkt','p,psia',2900,'supar', epsilon_up(i),'outp','massf','transport','mks','end');
            cea_data_iter = array2table([ceam_out.output.froz.aeat(3:length(ceam_out.output.froz.aeat)), ceam_out.output.froz.pressure(3:length(ceam_out.output.froz.pressure)), ...
                        ceam_out.output.froz.density(3:length(ceam_out.output.froz.density)), ceam_out.output.froz.temperature(3:length(ceam_out.output.froz.temperature)), ceam_out.output.froz.mach(3:length(ceam_out.output.froz.mach)), ...
                        ceam_out.output.froz.sonvel(3:length(ceam_out.output.froz.sonvel)), ceam_out.output.froz.isp(3:length(ceam_out.output.froz.isp)), ...
                        ceam_out.output.froz.cf(3:length(ceam_out.output.froz.cf)), ceam_out.output.froz.isp_vac(3:length(ceam_out.output.froz.isp_vac))]);
            cea_data_iter.Properties.VariableNames = ["aeat", "p", "rho", "T", "mach", "son", "isp", "cf", "ivac"];
            cea_data = [cea_data; cea_data_iter]; 
        end
        cea_data.p = cea_data.p*1e5; % Convert to Pascals

        %Plot our temperature and pressure as a function of x (finally)
        x_mod = xq(length(s)-(length(epsilon_down)+length(epsilon_up))+1:length(s));
        figure();
        plot(x_mod, cea_data.p/6895);
        xlabel("X (in)");
        ylabel("P (psia)");
        title("Fluid Pressure (Axially) in Solid Engine");
        grid on;

        figure();
        plot(x_mod, cea_data.T*1.8);
        xlabel("X (in)");
        ylabel("T (R)");
        title("Fluid Temperature (Axially) in Solid Engine");
        grid on;

        %Do T, P, rho all on the same plot relative to the chamber
        %conditions
        figure();
        plot(x_mod, cea_data.rho/cea_data.rho(1));
        hold on;
        plot(x_mod, cea_data.p/cea_data.p(1));
        hold on;
        plot(x_mod, cea_data.T/cea_data.T(1));
        xlabel("X (in)");
        ylabel("Stagnation Ratios");
        title("Stagnation Ratios Axially Through Solid Engine");
        legend(["Rho/Rho_o", "P/P_o", "T/T_o"]);
elseif(mode == "RDRE")
    %Engine Nominal Operation and Dimensions
    Weight = 350*9.81; % (N) [Upper Stage]
    Pc = 1.52e6; % (Pa)
    Dob = 4.0/39.37;  % (in -> m)
    Dib = 3.77/39.37; % (in -> m)
    Rt = (Dob-Dib)/2; % (m)
    Astar = (pi/4)*(Dob^2-Dib^2)
    R_throat_trad = sqrt(Astar/pi)
    R_curve = 0.382*R_throat_trad;
    optimize = false;
    
    if optimize == true
        %Run CEA with the given fuel and oxidizer for an RDRE engine and parse
        %ceam_out = CEA('reactants','fuel','N2H4','wt%',100,'t(k)',750,'oxid','N2O4(L)','wt%',100,'t(k)',298.15,'prob','rocket','case','RDRE','o/f', 1.4356, 'p,psia',300,'supar',25, 50, 75, 100,125,150,175,200,225,250,275,300,325,350,375,400,'equilibrium','output','mks','end','screen');
        fuel_percents = [0.6, 0.05, 0.3, 0.05]; % [HAN, AN, Methanol, Water] wt %
        fo_ratio = 50:95; % Fuel/Ox in 100%
        of_ratio = 1.1:0.025:2;
        Isp = [];
        
        for i = 1:length(of_ratio)
            ceam_out = CEA('reactants','fuel','N2H4','wt%',100,'t(k)',750,'oxid','N2O4(L)','wt%',100,'t(k)',298.15,'prob','rocket','case','RDRE','o/f', of_ratio(i), 'p,psia',300,'supar',25, 50, 75, 100,125,150,175,200,225,250,275,300,325,350,375,400,'equilibrium','output','mks','end','screen');
            ceam_out = CEA('reactants','fuel','N2H4','wt%',100,'t(k)',750,'oxid','N2O4(L)','wt%',100,'t(k)',298.15,'prob','rocket','case','RDRE','o/f', of_ratio(i), 'p,psia',300,'supar',25, 35, 45, 55, 65, 75, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150,'equilibrium','output','mks','end','screen');
            %ceam_out = CEA('reactants','name','HAN','wt',fuel_percents(1)*fo_ratio(i),'H',4,'N',2,'O',4,'h,kcal/mol',-95.300,'rho,g/cc',1.090,'name','AN','wt',fuel_percents(2)*fo_ratio(i),'H',4,'N',2,'O',3,'h,kcal/mol',-87.380,'rho,g/cc',1.725,'name','Methanol','wt',fuel_percents(3)*fo_ratio(i),'H',4,'C',1,'O',1,'h,kcal/mol',-57.0,'rho,g/cc',0.7918,'name','Water','wt',fuel_percents(4)*fo_ratio(i),'H',2,'O',1,'h,kcal/mol',-68.313,'rho,g/cc',1.000, 'name','N2O4(L)','wt%',100-fo_ratio(i),'t(k)',298.15, 'prob','rocket','case','RDRE', 'p,psia',300,'supar',25, 50, 75, 100,125,150,175,200,225,250,275,300,325,350,375,400,'equilibrium','output','mks','end','screen');
            %the presented results
            cea_data = array2table(squeeze([ceam_out.output.eql.aeat(3:length(ceam_out.output.eql.aeat)), ceam_out.output.eql.pressure(3:length(ceam_out.output.eql.pressure)), ...
                                ceam_out.output.eql.density(3:length(ceam_out.output.eql.density)), ceam_out.output.eql.mach(3:length(ceam_out.output.eql.mach)), ...
                                ceam_out.output.eql.sonvel(3:length(ceam_out.output.eql.sonvel)), ceam_out.output.eql.isp(3:length(ceam_out.output.eql.isp)), ...
                                ceam_out.output.eql.cf(3:length(ceam_out.output.eql.cf)), ceam_out.output.eql.isp_vac(3:length(ceam_out.output.eql.isp_vac))]));
            cea_data.Properties.VariableNames = ["aeat", "p", "rho", "mach", "son", "isp", "cf", "ivac"];
            cea_data.p = cea_data.p*1e5; % Convert to Pascals
            disp(cea_data.ivac);
    
            %Calculate mass flow and thrust
            mdot = Pc*Astar./(ceam_out.output.eql.cstar');
            Thrust = cea_data.cf*Pc*Astar; % (N)
        
            %Solve for corresponding RDRE bell nozzle lengths (considering plug)
            syms psi positive
            xi = zeros([length(cea_data.aeat), 1]);
            for i=1:length(xi)
                xi(i) = max(double(vpasolve(((Dib+2*psi)^2)/((Dob^2 - Dib^2))-cea_data.aeat(i) == 0, psi)));
            end
            Ln = xi./tand(15);
            length_trad = (R_throat_trad*(sqrt(cea_data.aeat)-1)+R_curve*((1/cosd(15))-1))/tand(15);
    
            %Set and store our Isp value for this iteration
            Isp = [Isp, cea_data.ivac];
        end
        %Plot our results for Isp v. Fuel Percentage 
        figure();
        surf(of_ratio, [25, 35, 45, 55, 65, 75, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150], Isp);
        xlabel("Fuel Percentage");
        ylabel("Expansion Ratio (Ae/At)")
        zlabel("Specific Impulse $Isp_{vac}$ (s)");
        title("Hydrazine Specific Impulse Performance");
    else
        ceam_out = CEA('reactants','fuel','N2H4','wt%',100,'t(k)',750,'oxid','N2O4(L)','wt%',100,'t(k)',298.15,'prob','rocket','case','RDRE','o/f', 1.325, 'p,psia',300,'supar',25, 50, 75, 100,125,150,175,200,225,250,275,300,325,350,375,400,'equilibrium','output','mks','end','screen');
        %ceam_out = CEA('reactants','fuel','N2H4','wt%',100,'t(k)',750,'oxid','N2O4(L)','wt%',100,'t(k)',298.15,'prob','rocket','case','RDRE','o/f', of_ratio(i), 'p,psia',300,'supar',25, 35, 45, 55, 65, 75, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150,'equilibrium','output','mks','end','screen');
        %ceam_out = CEA('reactants','name','HAN','wt',fuel_percents(1)*fo_ratio(i),'H',4,'N',2,'O',4,'h,kcal/mol',-95.300,'rho,g/cc',1.090,'name','AN','wt',fuel_percents(2)*fo_ratio(i),'H',4,'N',2,'O',3,'h,kcal/mol',-87.380,'rho,g/cc',1.725,'name','Methanol','wt',fuel_percents(3)*fo_ratio(i),'H',4,'C',1,'O',1,'h,kcal/mol',-57.0,'rho,g/cc',0.7918,'name','Water','wt',fuel_percents(4)*fo_ratio(i),'H',2,'O',1,'h,kcal/mol',-68.313,'rho,g/cc',1.000, 'name','N2O4(L)','wt%',100-fo_ratio(i),'t(k)',298.15, 'prob','rocket','case','RDRE', 'p,psia',300,'supar',25, 50, 75, 100,125,150,175,200,225,250,275,300,325,350,375,400,'equilibrium','output','mks','end','screen');
        %the presented results
        cea_data = array2table(squeeze([ceam_out.output.eql.aeat(3:length(ceam_out.output.eql.aeat)), ceam_out.output.eql.pressure(3:length(ceam_out.output.eql.pressure)), ...
                            ceam_out.output.eql.density(3:length(ceam_out.output.eql.density)), ceam_out.output.eql.mach(3:length(ceam_out.output.eql.mach)), ...
                            ceam_out.output.eql.sonvel(3:length(ceam_out.output.eql.sonvel)), ceam_out.output.eql.isp(3:length(ceam_out.output.eql.isp)), ...
                            ceam_out.output.eql.cf(3:length(ceam_out.output.eql.cf)), ceam_out.output.eql.isp_vac(3:length(ceam_out.output.eql.isp_vac))]));
        cea_data.Properties.VariableNames = ["aeat", "p", "rho", "mach", "son", "isp", "cf", "ivac"];
        cea_data.p = cea_data.p*1e5; % Convert to Pascals
        disp(cea_data.ivac);

        %Calculate mass flow and thrust
        mdot = Pc*Astar./(ceam_out.output.eql.cstar');
        Thrust = cea_data.ivac*mdot*9.81; % (N)
    
        %Solve for corresponding RDRE bell nozzle lengths (considering plug)
        syms psi positive
        xi = zeros([length(cea_data.aeat), 1]);
        for i=1:length(xi)
            xi(i) = max(double(vpasolve(((Dib+2*psi)^2)/((Dob^2 - Dib^2))-cea_data.aeat(i) == 0, psi)));
        end
        Ln = xi./tand(15);
        length_trad = (R_throat_trad*(sqrt(cea_data.aeat)-1)+R_curve*((1/cosd(15))-1))/tand(15);
        
    end

    %Plot our results for length distributions
    figure();
    plot(cea_data.aeat, Ln*3.28084);
    hold on;
    plot(cea_data.aeat, length_trad*3.28084);
    xlabel("Ae/At");
    ylabel("Nozzle Length $L_n$ (ft)");
    title("Nozzle Length v. Area Ratio");
    legend(["RDRE", "Traditional"]);
     
    %Plot our results for Isp vs area ratio
    figure();
    plot(cea_data.aeat, cea_data.ivac);
    xlabel("Ae/At");
    ylabel("Specific Impulse $Isp_{vac}$ (s)");
    title("Specific Impulse (Vacuum) v. Area Ratio");
    
    %Plot results for Thrust
    figure();
    plot(cea_data.aeat, Thrust/Weight);
    xlabel("Ae/At");
    ylabel("Thrust/Weight");
    title("Thrust to Weight v. Area Ratio");

    %Plot our results for Engine Radius
    figure();
    plot(cea_data.aeat, sqrt(cea_data.aeat*Astar/pi)*39.37);
    xlabel("Ae/At");
    ylabel("Engine Radius (in)");
    title("Engine Radius v. Area Ratio");
end

