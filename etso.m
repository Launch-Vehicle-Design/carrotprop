%% Engine and Tank System Optimization- Nathan Rand
% 01/29/2024
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

rocket_tank_pressure_simulation;

function rocket_tank_pressure_simulation
    % Load in trajectory
    load('alti.mat');
    
    % Constants and Initial Conditions
    P0 = 2e6; % Initial pressure in the tank at ground level (Pa)
    V0 = 0.105; % Volume of the tank (m^3)
    T0 = 293; % Initial temperature of the propellant (K), assuming sea level standard temperature

    % Time steps for each phase (in seconds)
    phase1_duration = 5 * 60; % F15E climb: 5 mins
    phase2_duration = 35 * 60; % F15E cruise: 45 mins
    phase3_duration = 250; % 1st stage climb: 4 mins and 10 s
    phase4_duration = 250; % 2nd stage climb and orbit: 4 mins and 10 s

    % Altitude function for each phase
    phase1_altitude = @(t) 12000 * t / phase1_duration; % Climbing to 12,000 m in 5 mins
    phase2_altitude = @(t) 12000; % Constant altitude of 12,000 m during cruise
    phase3_altitude = alti; % Climbing to 400,000 ft (121920 m) in 10 mins
    phase4_altitude = alti; % Climbing to 400,000 ft (121920 m) in 10 mins

    % Total simulation time and time steps
    total_time = phase1_duration + phase2_duration + phase3_duration + phase4_duration;
    time_steps = 0:1:total_time; % Time steps every 10 seconds

    % Initialize arrays for temperature and pressure
    T = zeros(size(time_steps)); % Propellant temperature
    T_ext = zeros(size(time_steps)); % Atmospheric temperature for comp.
    P = zeros(size(time_steps)); % Tank pressure
    Power = zeros(size(time_steps));

    % Loop over time steps to calculate temperature and pressure
    for i = 1:length(time_steps)
        t = time_steps(i);
        
        % Determine current altitude based on phase
        if t <= phase1_duration
            altitude = phase1_altitude(t);
        elseif t <= (phase1_duration + phase2_duration)
            altitude = phase2_altitude(t - phase1_duration);
        elseif t <= (phase1_duration + phase2_duration + phase3_duration)
            altitude = phase3_altitude(t-(phase2_duration+phase1_duration));
        else
            altitude = phase4_altitude(t-(phase2_duration+phase1_duration));
        end

        % Check if altitude is within the range of atmosisa
        T_ext(i) = get_external_temperature(altitude);

        % Simple heat transfer model to update propellant temperature
        T_o = T_ext(i);
        if i > 1 % After the first iteration, assign our initial prop temp to the previous iter.
            T_o = T(i-1);
        end
        [T(i), Power(i)] = heat_transfer_model(T_ext(i), T_o, altitude, t);

        % Calculate pressure in the tank using ideal gas law
        P(i) = P0 * (T(i)/T0);
    end

    % Plot results
    figure;
    plot(time_steps/60, P*145/1e6);
    xlabel('Time (minutes)');
    ylabel('Pressure (psia)');
    title('Tank Pressure Over Time');
    
    figure;
    plot(time_steps/60, T*1.8);
    hold on;
    yline(495, LineWidth=2, LineStyle="--"); % (K) [Freezing Point of Hydrazine]
    xlabel('Time (minutes)');
    ylabel('Temperature (R)');
    title('Propellant Temperature Over Time');
    legend(["Temp. With Heater", "Freezing Point"]);
    ylim([495*0.98, max(T)*1.8*1.02]);

    figure;
    plot(time_steps/60, Power*3.41);
    xlabel('Time (minutes)');
    ylabel('Heater Required Power (BTU)');
    title('Heater Required Power Over Time');
end

function [T_prop, P_heater] = heat_transfer_model(T_external, T_o, altitude, time_step)
    % Constants
    lower_transition_altitude = 65000; % Transition begins at 65 km
    upper_transition_altitude = 85000; % Transition ends at 85 km
    Cp = 2530; % (J/kgK)
    m = 105.06; % (kg)

    % Linear interpolation for transition weights
    if altitude > lower_transition_altitude && altitude < upper_transition_altitude
        transition_weight = (altitude - lower_transition_altitude) / (upper_transition_altitude - lower_transition_altitude);
    elseif altitude >= upper_transition_altitude
        transition_weight = 1;
    else
        transition_weight = 0;
    end

    % Convective/conductive heat transfer component
    Q_convective = conductive_convective_heat_transfer(T_external, T_o, time_step) * (1 - transition_weight);
    
    % Radiative heat transfer component
    Q_radiative = radiative_heat_transfer(T_external, T_o, time_step) * transition_weight;
    
    % Total heat transfer is the sum of convective and radiative components
    Q_total = Q_convective + Q_radiative;

    % Calculate the initial temperature change without heater
    deltaT = Q_total * time_step / (Cp * m);
    T_without_heater = T_o + deltaT;
    T_wo_heater = T_without_heater;

    % Check if external heating is needed
    freezing_temp_margin = 285.0; % (K)
    if T_without_heater < freezing_temp_margin
        % Calculate additional heat required to reach the threshold
        Q_needed = (freezing_temp_margin - T_without_heater) * Cp * m / time_step;
        
        % Assuming the heater's power input is directly converted to heat
        P_heater = Q_needed / time_step; % Power in Watts
        
        % Update the temperature with the heater's contribution
        T_prop = freezing_temp_margin;
    else
        P_heater = 0; % No additional heating required
        T_prop = T_without_heater;
    end
end

function Q_conv = conductive_convective_heat_transfer(T_external, T_o, time_step)
    % Define thermal conductivities (in W/mK)
    k_air = 0.024;
    k_carbon_fiber = 100;
    k_aluminum = 210;
    k_N2H4 = 0.32;

    % Define heat transfer coefficients (in W/m^2K)
    h_air = 10; % Convective heat transfer coefficient for air
    h_N2H4 = 100; % Convective heat transfer coefficient for N2H4

    % Thicknesses of the materials (in meters)
    t_carbon_fiber = 1e-3;
    t_aluminum = 1e-3;
    
    %Surface area of the tank (in m^2)
    A = 1.3; % For fuel tank

    % Calculate Thermal Resistances
    R_conv_air = 1 / (h_air * A); % Convective resistance for air
    R_cond_carbon = t_carbon_fiber / (k_carbon_fiber * A); % Conduction resistance for carbon fiber
    R_cond_aluminum = t_aluminum / (k_aluminum * A); % Conduction resistance for aluminum
    R_conv_N2H4 = 1 / (h_N2H4 * A); % Convective resistance for N2H4

    % Total Resistance
    R_total = R_conv_air + R_cond_carbon + R_cond_aluminum + R_conv_N2H4;

    % Heat Transfer (Q = deltaT / R_total)
    deltaT = T_external - T_o;
    Q_conv = deltaT / R_total;
end

function Q_rad = radiative_heat_transfer(T_external, T_o, time_step)
    sigma = 5.67e-8; % Stefan-Boltzmann constant in W/m^2K^4
    solar_constant = 1361; % Solar constant in W/m^2
    A = 1.3; % m^2
    absorptivity = 0.2; % Ability to absorb solar radiation
    epsilon = 0.85;

    % Calculate absorbed solar radiation
    Q_absorbed = solar_constant * absorptivity * A;

    % Calculate radiative cooling
    Q_radiated = epsilon * sigma * A * T_o^4;

    % Net heat transfer to the tank
    Q_rad = Q_absorbed - Q_radiated;
end

function T_external = get_external_temperature(altitude)
    if altitude <= 85000 % Use atmosisa for altitudes up to 85 km
        [T_external, ~, ~, ~] = atmosisa(altitude);
    else % Use a simple model for the thermosphere
        T_external = thermosphere_temperature_model(altitude);
    end
end

function T_ext = thermosphere_temperature_model(altitude)
    % Define temperature ranges and corresponding altitudes in the thermosphere
    base_altitude = 85000; % 85 km in meters
    max_altitude = 500000; % 500 km in meters

    % Define temperatures at base and max altitudes (in Kelvin)
    base_temp = 186.946; % (K)
    max_temp = 2273; % (K)

    if altitude <= base_altitude
        T_ext = base_temp;
    elseif altitude >= max_altitude
        T_ext = max_temp;
    else
        % Linear interpolation between base_temp and max_temp
        T_ext = base_temp + ((max_temp - base_temp) / (max_altitude - base_altitude)) * (altitude - base_altitude);
    end
end