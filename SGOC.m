%% Solid Grain Optimization Code - Nathan Rand
% 02/03/2024
close all;
clc;
clear;

% Global variables to store optimization history (bad I know but meh)
global optimHistory;
optimHistory.x = [];
optimHistory.isp = [];
optimHistory.mode = "COTS";
mode = optimHistory.mode;
optimHistory.brs = "None";
brs = optimHistory.brs;
if mode == "CL20"
    % Equality constraint for cl20 + gap + ap = 100
    Aeq = [1, 1, 1, 1, 1];
    beq = 100;
    
    % Bounds
    lb = [17.5, 12, 0, 0, 10]; % Lower bounds for cl20, gap, and ap
    ub = [100, 20, 100, 20, 10]; % Upper bounds for cl20, gap, and ap
    
    % Initial guess
    x0 = [42.5, 15, 27.5, 5, 12.5]; % Example: 42.5% cl20, 15% gap, 42.5% ap
    
    % Set options for fmincon
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', 'MaxFunctionEvaluations', 1e4, 'MaxIterations', 1e4, 'StepTolerance', 1e-6, 'OptimalityTolerance', 1e-6, 'OutputFcn', @outfun);
    
    % Run the optimization
    [x_opt, fval_opt] = fmincon(@objectiveFunction, x0, [], [], Aeq, beq, lb, ub, [], options);
    disp(x_opt);
    disp(fval_opt);
    
    % Plotting results after optimization
    figure();
    subplot(2,1,1);
    plot(optimHistory.isp);
    xlabel('Iteration');
    ylabel('Isp (s)');
    title('Isp Optimization Progress');
    
    subplot(2,1,2);
    plot(optimHistory.x);
    xlabel('Iteration');
    ylabel('Component Concentration (\%)');
    legend('CL20', 'GAP', 'AP', 'AL', "Oxamide");
    title('Propellant Component Concentrations Progress');
elseif mode == "COTS"
    if brs == "Oxamide"
        % Equality constraint for cl20 + gap + ap = 100
        Aeq = [1, 1, 1, 1, 1];
        beq = 100;
        
        % Bounds
        lb = [12, 0, 0, 10, 5]; 
        ub = [20, 100, 20, 10, 20]; 
        
        % Initial guess
        x0 = [15, 50, 10, 10, 5]; % Example: HTPB, AP, AL, Oxamide, HMX
    elseif brs == "None"
        % Equality constraint for cl20 + gap + ap = 100
        Aeq = [1, 1, 1];
        beq = 100;
        
        % Bounds
        lb = [12, 0, 0]; 
        ub = [20, 100, 20]; 
        
        % Initial guess
        x0 = [15, 50, 10]; % Example: HTPB, AP, AL
    end
    
    % Set options for fmincon
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', 'MaxFunctionEvaluations', 1e4, 'MaxIterations', 1e4, 'StepTolerance', 1e-6, 'OptimalityTolerance', 1e-6, 'OutputFcn', @outfun);
    
    % Run the optimization
    if mode == "CL20"
        [x_opt, fval_opt] = fmincon(@objectiveFunction, x0, [], [], Aeq, beq, lb, ub, [], options);
    elseif mode == "COTS"
        if brs == "Oxamide"
            [x_opt, fval_opt] = fmincon(@objectiveFunctionCOTS_brs, x0, [], [], Aeq, beq, lb, ub, [], options);
        elseif brs == "None"
            [x_opt, fval_opt] = fmincon(@objectiveFunctionCOTS, x0, [], [], Aeq, beq, lb, ub, [], options);
        end
    end
    disp(x_opt);
    disp(fval_opt);
    
    % Plotting results after optimization
    figure();
    subplot(2,1,1);
    plot(optimHistory.isp);
    xlabel('Iteration');
    ylabel('Isp (s)');
    title('Isp Optimization Progress');
    
    subplot(2,1,2);
    plot(optimHistory.x);
    xlabel('Iteration');
    ylabel('Component Concentration');
    if mode == "CL20"
        legend('CL20', 'GAP', 'AP', 'AL', "Oxamide");
    elseif mode == "COTS"
        if brs == "Oxamide"
            legend('HTPB', 'AP', 'AL', "Oxamide", "HMX");
        elseif brs == "None"
            legend('HTPB', 'AP', 'AL');
        end
    end
    title('Propellant Component Concentrations Progress');

end

function stop = outfun(x,optimValues,state)
    global optimHistory;
    stop = false;
    switch state
        case 'iter'
            optimHistory.x = [optimHistory.x; x];
            optimHistory.isp = [optimHistory.isp; -optimValues.fval];
    end
end

function negIspEff = objectiveFunction(x)
    % Unpack the input array
    cl20 = x(1);
    gap = x(2);
    ap = x(3);
    al = x(4);
    oxamide = x(5);

    %Define problem wide variables
    H_release = 12192; % (m)
    H_orbit = 251.46e3; % (m)
    Weight = 2267.962*9.81; % (N)
    g = 9.81; % m/s^2
    Ru = 8314; % (J/kgK)
    [Ta, a, Pa, rho] = atmosisa(H_release);

    %Process Results through CEA
    try
        ceam_out = CEA('reac','name','CL20','C', 6.28, 'H', 6.20, 'N', 12.55, 'O', 12, 'h,Kj/mol', 581.87,'wt%', cl20, 't(k)',298.15, 'name','GAP','C', 388, 'H', 330, 'O', 4, 'N', 192, 'h,Kj/mol', 6060.0, 'wt%', gap, 't(k)', 298.15, 'name', 'NH4CLO4(I)', 't(k)', 298.15, 'wt(%)', ap, 'name', 'AL', 't(k)', 298.15, 'wt(%)', al, 'name', 'Oxamide', 'C', 2, 'H', 4, 'O', 2, 'N', 2, 'h,Kj/mol', -364.8, 't(k)', 298.15, 'wt(%)', oxamide, 'prob','rkt','p,psia',1088,'supar',5,10,15,20,25,30,35,40,45,50,55,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,300,'outp','massf','transport','mks','end');
        cea_data = array2table(squeeze([ceam_out.output.froz.aeat(3:length(ceam_out.output.froz.aeat)), ceam_out.output.froz.pressure(3:length(ceam_out.output.froz.pressure)), ...
                            ceam_out.output.froz.density(3:length(ceam_out.output.froz.density)), ceam_out.output.froz.temperature(3:length(ceam_out.output.froz.temperature)), ceam_out.output.froz.mach(3:length(ceam_out.output.froz.mach)), ...
                            ceam_out.output.froz.sonvel(3:length(ceam_out.output.froz.sonvel)), ceam_out.output.froz.isp(3:length(ceam_out.output.froz.isp)), ...
                            ceam_out.output.froz.cf(3:length(ceam_out.output.froz.cf)), ceam_out.output.froz.isp_vac(3:length(ceam_out.output.froz.isp_vac))])');
        cea_data.Properties.VariableNames = ["aeat", "p", "rho", "T", "mach", "son", "isp", "cf", "ivac"];
        cea_data.p = cea_data.p*1e5; % Convert to Pascals
    catch
        negIspEff = 999; 
        return
    end
    %Chamber Pressure: 7.5 MPa
    Pc = 7.5e6; % (Pa)
    Rt = 0.07; % meters (approx. 0.72 inches)
    rho = 1948.8; % (kg/m^3)
    rb = 23.5/1e3; % (m/s)
    Astar = pi*Rt^2; % Throat Area (m^2)
    R_curve = 0.382*Rt; % meters
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
    disp(length_trunc(find(Isp_percent==0.96)));
    disp(epsilon(find(Isp_percent==0.96)));
    
    % Return the negative of effective Isp for minimization
    negIspEff = -0.96*max(Isp_eff);
end

function negIspEff = objectiveFunctionCOTS(x)
    % Unpack the input array
    htpb = x(1);
    ap = x(2);
    al = x(3);

    %Define problem wide variables
    H_release = 12192; % (m)
    H_orbit = 251.46e3; % (m)
    Weight = 2267.962*9.81; % (N)
    g = 9.81; % m/s^2
    Ru = 8314; % (J/kgK)
    [Ta, a, Pa, rho] = atmosisa(H_release);
    %Chamber Pressure: 9.5 MPa
    Pc = 9.5; % (Pa)

    %Process Results through CEA
    try
        ceam_out = CEA('reac','name','C4H6,butadiene' ,'wt(%)', htpb, 't(k)', 298.15, 'name', 'NH4CLO4(I)', 't(k)', 298.15, 'wt(%)', ap, 'name', 'AL(cr)', 't(k)', 298.15, 'wt(%)', al, 'prob','rkt','p,psia',Pc*145.038,'supar',5,10,15,20,25,30,35,40,45,50,55,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,300,'outp','massf','transport','mks','end');
        cea_data = array2table(squeeze([ceam_out.output.froz.aeat(3:length(ceam_out.output.froz.aeat)), ceam_out.output.froz.pressure(3:length(ceam_out.output.froz.pressure)), ...
                            ceam_out.output.froz.density(3:length(ceam_out.output.froz.density)), ceam_out.output.froz.temperature(3:length(ceam_out.output.froz.temperature)), ceam_out.output.froz.mach(3:length(ceam_out.output.froz.mach)), ...
                            ceam_out.output.froz.sonvel(3:length(ceam_out.output.froz.sonvel)), ceam_out.output.froz.isp(3:length(ceam_out.output.froz.isp)), ...
                            ceam_out.output.froz.cf(3:length(ceam_out.output.froz.cf)), ceam_out.output.froz.isp_vac(3:length(ceam_out.output.froz.isp_vac))])');
        cea_data.Properties.VariableNames = ["aeat", "p", "rho", "T", "mach", "son", "isp", "cf", "ivac"];
        cea_data.p = cea_data.p*1e5; % Convert to Pascals
    catch
        negIspEff = 999; 
        return
    end
    Rt = 0.047/2; % meters (approx. 0.72 inches)
    rho = 1948.8; % (kg/m^3)
    rb = 23.5/1e3; % (m/s)
    Astar = pi*Rt^2; % Throat Area (m^2)
    R_curve = 0.382*Rt; % meters
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
    disp(length_trunc(find(Isp_percent==0.9525))*3.281);
    disp(epsilon(find(Isp_percent==0.9525)));
    
    % Return the negative of effective Isp for minimization
    negIspEff = -0.9525*max(Isp_eff);
end


function negIspEff = objectiveFunctionCOTS_brs(x)
    % Unpack the input array
    htpb = x(1);
    ap = x(2);
    al = x(3);
    oxamide = x(4);
    hmx = x(5);

    %Define problem wide variables
    H_release = 12192; % (m)
    H_orbit = 251.46e3; % (m)
    Weight = 2267.962*9.81; % (N)
    g = 9.81; % m/s^2
    Ru = 8314; % (J/kgK)
    [Ta, a, Pa, rho] = atmosisa(H_release);

    %Process Results through CEA
    try
        ceam_out = CEA('reac','name','HMX','C', 4, 'H', 8, 'N', 8, 'O', 8, 'h,Kj/mol', 75.0,'wt%', hmx, 't(k)',298.15, 'name','C4H6,butadiene' ,'wt(%)', htpb, 't(k)', 298.15, 'name', 'NH4CLO4(I)', 't(k)', 298.15, 'wt(%)', ap, 'name', 'AL(cr)', 't(k)', 298.15, 'wt(%)', al, 'name', 'Oxamide', 'C', 2, 'H', 4, 'O', 2, 'N', 2, 'h,Kj/mol', -364.8, 't(k)', 298.15, 'wt(%)', oxamide, 'prob','rkt','p,psia',1088,'supar',5,10,15,20,25,30,35,40,45,50,55,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,300,'outp','massf','transport','mks','end');
        cea_data = array2table(squeeze([ceam_out.output.froz.aeat(3:length(ceam_out.output.froz.aeat)), ceam_out.output.froz.pressure(3:length(ceam_out.output.froz.pressure)), ...
                            ceam_out.output.froz.density(3:length(ceam_out.output.froz.density)), ceam_out.output.froz.temperature(3:length(ceam_out.output.froz.temperature)), ceam_out.output.froz.mach(3:length(ceam_out.output.froz.mach)), ...
                            ceam_out.output.froz.sonvel(3:length(ceam_out.output.froz.sonvel)), ceam_out.output.froz.isp(3:length(ceam_out.output.froz.isp)), ...
                            ceam_out.output.froz.cf(3:length(ceam_out.output.froz.cf)), ceam_out.output.froz.isp_vac(3:length(ceam_out.output.froz.isp_vac))])');
        cea_data.Properties.VariableNames = ["aeat", "p", "rho", "T", "mach", "son", "isp", "cf", "ivac"];
        cea_data.p = cea_data.p*1e5; % Convert to Pascals
    catch
        negIspEff = 999; 
        return
    end
    %Chamber Pressure: 7.5 MPa
    Pc = 7.5e6; % (Pa)
    Rt = 0.05; % meters (approx. 0.72 inches)
    rho = 1948.8; % (kg/m^3)
    rb = 23.5/1e3; % (m/s)
    Astar = pi*Rt^2; % Throat Area (m^2)
    R_curve = 0.382*Rt; % meters
    %Determine area ratios that are optimal for specific Pa values
    Isp_40 = (Astar*(cea_data.rho.*cea_data.aeat.*(cea_data.mach.*cea_data.son).^2 + cea_data.aeat.*(cea_data.p-Pa)))./(Astar.*cea_data.rho.*cea_data.aeat.*(cea_data.mach.*cea_data.son).*g);
    Isp_eff = 0.98*(Isp_40 + (2/3)*(cea_data.ivac - Isp_40));
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
    disp(length_trunc(find(Isp_percent==0.96)));
    disp(epsilon(find(Isp_percent==0.96)));
    
    % Return the negative of effective Isp for minimization
    negIspEff = -0.96*max(Isp_eff);
end



