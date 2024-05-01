function [prop,nozz,var,perf] = internalballistics_V1(param,CEA,prop,geom,nozz,ignit,flag)
%INTERNALBALLISTICS Summary of this function goes here
% INPUTS
%   "param" includes the following: temps [K], sig_p [/K], h [m], dt [s]
%   "CEA"   includes the following: tc [K], cstar [m/s], gamma, mw [kg/kmol]
%   "prop   includes the following: yi [%], rho_i [kg/m^3], rb [m/s], pc [MPa]
%   "geom"  includes the following: ID [m], OD [m], IL [m], type ["string"]
%   "nozz"  includes the following: at [m^2], ae [m^2], edot [m/s], angle [deg]
%   "ignit" includes the following: pc [Pa], tc [K], vol_c [m^3]
%   "flag"  includes the following: nozz_erosion, flag1, flag2
% OUTPUTS

% define additional constants
param.g0 = 9.81; % m/s^2 -> can add changing gravity functionality
param.Ru = 8.31446261815324e3; % [J/(kmol*K)]

% calculate R
CEA.R = param.Ru/CEA.mw;

% calculate mixture density
prop.rho = 1./sum(prop.yi./prop.rho); % mixture density

% calculate initial propellant mass
if geom.type == "end-burner"
    prop.mp = (pi*geom.OD^2/4)*geom.IL*prop.rho; % [kg]
elseif geom.type == "center-perf"
    prop.mp = (pi*(geom.OD^2/4-geom.ID^2/4))*geom.IL*prop.rho;
elseif geom.type == "rod & tube"
    prop.mp = (pi/4*(geom.OD^2+geom.D1^2-geom.D2^2))*geom.IL*prop.rho;
elseif geom.type == "bates"
    prop.mp = (pi*(geom.OD^2/4-geom.ID^2/4))*geom.IL*prop.rho;
end

% solve for a and n
prop.pc_log = log(prop.pc);
prop.rb_log = log(prop.rb);
p = polyfit(prop.pc_log,prop.rb_log,1);
prop.a = exp(p(2));
prop.n = p(1);

% calculate additional nozzle parameters
nozz.eps = nozz.ae/nozz.at;
nozz.dt = sqrt(4*nozz.at/pi);
nozz.de = sqrt(4*nozz.ae/pi);

% evaluate atmospheric conditions
[atmo.ta, atmo.a, atmo.pa, atmo.rho] = atmosisa(param.h);

% initialize variables
var.pc = zeros(2e5,1);
var.tc = zeros(2e5,1);
var.vol_c = zeros(2e5,1);
if flag.nozz_erosion
    var.dt = zeros(2e5,1);
    var.at = zeros(2e5,1);
end
var.rb = zeros(2e5,1);
var.mdot_gen = zeros(2e5,1);
var.mdot_out = zeros(2e5,1);
% var.mdot_sto = zeros(2e5,1);
var.dv_dt = zeros(2e5,1);
var.dpc_dt = zeros(2e5,1);
var.mp_burnt = zeros(2e5,1);
var.length = zeros(2e5,1);
var.time = zeros(2e5,1);

% initialize performance parameters
perf.me = zeros(2e5,1);
perf.pe = zeros(2e5,1);
perf.cf = zeros(2e5,1);
perf.f  = zeros(2e5,1);

% initialize starting conditions
var.pc(1) = ignit.pc; % Pa  -> pressure after ignition
var.tc(1) = ignit.tc; % K   -> temperature after ignition
var.vol_c(1) = ignit.vol_c; % m^3 -> initial volume of combustion chamber to throat

% calculate initial burn area
if geom.type == "center-perf"
    var.ab = zeros(2e5,1);
    var.ab(1) = 2*pi*(geom.ID/2)*geom.IL;
    var.id = zeros(2e5,1);
    var.id(1) = geom.ID;
end
if geom.type == "end-burner" && flag.insulation
    var.texp = zeros(2e5,1);
    var.texp(1) = param.tb;
    var.th_ins = zeros(2e5,1);
    var.th_ins(1) = var.texp(1)*param.edot*param.fs;
    var.id = zeros(2e5,1);
    var.id(1) = geom.OD - 2*var.th_ins(1);
    var.ab = zeros(2e5,1);
    var.ab(1) = pi*(var.id(1)^2/4);
elseif geom.type == "end-burner"
    var.id = zeros(2e5,1);
    var.id(1) = geom.OD;
    var.ab = zeros(2e5,1);
    var.ab(1) = pi*(var.id(1)^2/4);
end
if geom.type == "end-burner" && flag.dome
    var.height_dome = zeros(2e5,1);
end
if geom.type == "rod & tube"
    var.D1 = zeros(2e5,1);
    var.D1(1) = geom.D1;
    var.D2 = zeros(2e5,1);
    var.D2(1) = geom.D2;
    var.ab = zeros(2e5,1);
    var.ab(1) = pi*geom.IL*(var.D1(1)+var.D2(1));
end
if geom.type == "bates"
    var.ab = zeros(2e5,1);
    var.ab(1) = 2*pi*(geom.ID/2)*geom.IL + 2*pi*(geom.OD^2/4 - geom.ID^2/4);
    var.id = zeros(2e5,1);
    var.id(1) = geom.ID;
end

if flag.nozz_erosion
    var.dt(1) = nozz.dt;
    var.at(1) = nozz.at;
    var.eps(1) = nozz.ae/nozz.at;
end

i = 2;
while var.mp_burnt(i-1) < 1278.8174
    % calculate burn rate
    if flag.temp_sens_high
        var.rb(i) = prop.a*(var.pc(i-1)*1e-6)^prop.n*exp(param.sig_p*(param.temps(3)-param.temps(2))); % calculate current burn rate
    end
    if flag.temp_sens_low
        var.rb(i) = prop.a*(var.pc(i-1)*1e-6)^prop.n*exp(param.sig_p*(param.temps(1)-param.temps(2))); % calculate current burn rate
    end
    if flag.temp_sens_ambient
        var.rb(i) = prop.a*(var.pc(i-1)*1e-6)^prop.n; % calculate current burn rate
    end

    % calculate mdot_gen
    var.mdot_gen(i) = var.ab(i-1)*var.rb(i)*prop.rho;
    
    % calculate mdot_out
    if flag.nozz_erosion
        var.mdot_out(i) = var.pc(i-1)*var.at(i-1)/CEA.cstar;
    else
        var.mdot_out(i) = var.pc(i-1)*nozz.at/CEA.cstar;
    end

    % var.mdot_sto(i) = var.mdot_gen(i) - var.mdot_out(i);
    var.dv_dt(i) = var.mdot_gen(i)/prop.rho;
    var.dpc_dt(i) = CEA.R*var.tc(i-1)/var.vol_c(i-1)*(var.mdot_gen(i) - var.mdot_out(i)) - (var.pc(i-1)/var.vol_c(i-1))*var.dv_dt(i);
    
    % calcualate burnt propellant mass and consumed grain dimensions
    var.mp_burnt(i) = var.mp_burnt(i-1) + var.mdot_gen(i)*param.dt; % kg
    if geom.type == "center-perf"
        var.id(i) = var.id(i-1) + 2*var.rb(i)*param.dt; % m
    end
    if geom.type == "end-burner" && flag.insulation
        var.length(i) = var.length(i-1) + var.rb(i)*param.dt; % m
        if var.time <= param.tb
            var.texp(i) = param.tb - var.time(i-1);
            var.th_ins(i) = var.texp(i)*param.edot*param.fs;
            var.id(i) = geom.OD - 2*var.th_ins(i);
        else
            var.texp(i) = var.texp(i-1);
            var.th_ins(i) = var.th_ins(i-1);
            var.id(i) = var.id(i-1);
        end
    elseif geom.type == "end-burner"
        var.id(i) = var.id(i-1);
        var.length(i) = var.length(i-1) + var.rb(i)*param.dt; % m
    end
    if geom.type == "end-burner" && flag.dome 
        var.height_dome(i) = 0; % m
    end
    if geom.type == "end-burner" && flag.dome && var.length(i) > geom.IL
        var.height_dome(i) = var.height_dome(i-1) + var.rb(i)*param.dt; % m
    end
    if geom.type == "rod & tube"
        var.D1(i) = var.D1(i-1) - 2*var.rb(i)*param.dt;
        var.D2(i) = var.D2(i-1) + 2*var.rb(i)*param.dt;
    end
    if geom.type == "bates"
        var.id(i) = var.id(i-1) + 2*var.rb(i)*param.dt; % m
        var.length(i) = var.length(i-1) + 2*var.rb(i)*param.dt; % m
    end 
    % calculate performance parameters
    func = @(me) 1/me*sqrt((2/(CEA.gamma+1)*(1+(CEA.gamma-1)/2*me^2))^((CEA.gamma+1)/(CEA.gamma-1)));
    if flag.nozz_erosion
        eq1 = @(me) func(me) - var.eps(i-1);
    else
        eq1 = @(me) func(me) - nozz.eps;
    end
    
    perf.me(i) = fzero(eq1,1);
    perf.pe(i) = var.pc(i-1)/(1+(CEA.gamma-1)/2*perf.me(i)^2)^(CEA.gamma/(CEA.gamma-1));
    
    if flag.nozz_eta
        if flag.nozz_erosion
            perf.cf(i) = nozz.lambda*sqrt(2*CEA.gamma^2/(CEA.gamma-1)*(2/(CEA.gamma+1))^((CEA.gamma+1)/(CEA.gamma-1))*(1-(perf.pe(i)/var.pc(i-1))^((CEA.gamma-1)/CEA.gamma))) + (perf.pe(i)-atmo.pa)/var.pc(i-1)*var.eps(i-1);
        else
            perf.cf(i) = nozz.lambda*sqrt(2*CEA.gamma^2/(CEA.gamma-1)*(2/(CEA.gamma+1))^((CEA.gamma+1)/(CEA.gamma-1))*(1-(perf.pe(i)/var.pc(i-1))^((CEA.gamma-1)/CEA.gamma))) + (perf.pe(i)-atmo.pa)/var.pc(i-1)*nozz.eps;
        end
    else
        if flag.nozz_erosion
            perf.cf(i) = sqrt(2*CEA.gamma^2/(CEA.gamma-1)*(2/(CEA.gamma+1))^((CEA.gamma+1)/(CEA.gamma-1))*(1-(perf.pe(i)/var.pc(i-1))^((CEA.gamma-1)/CEA.gamma))) + (perf.pe(i)-atmo.pa)/var.pc(i-1)*var.eps(i-1);
        else
            perf.cf(i) = sqrt(2*CEA.gamma^2/(CEA.gamma-1)*(2/(CEA.gamma+1))^((CEA.gamma+1)/(CEA.gamma-1))*(1-(perf.pe(i)/var.pc(i-1))^((CEA.gamma-1)/CEA.gamma))) + (perf.pe(i)-atmo.pa)/var.pc(i-1)*nozz.eps;
        end
    end

    if flag.nozz_erosion
        perf.f(i)  = perf.cf(i)*var.pc(i-1)*var.at(i-1);
    else
        perf.f(i)  = perf.cf(i)*var.pc(i-1)*nozz.at;
    end
    
    % calculate pressure increase
    var.pc(i) = var.pc(i-1) + var.dpc_dt(i)*param.dt;

    % calculate conditions for next interation
    var.tc(i) = CEA.tc; % need to integrate CEA for the current chamber pressure
    
    % calculate volume of combustion chamber
    var.vol_c(i) = var.vol_c(i-1) + var.dv_dt(i)*param.dt;
    
    % calculate new burn area
    if geom.type == "center-perf"
        var.ab(i) = 2*pi*(var.id(i)/2)*geom.IL;
    elseif geom.type == "end-burner"
        var.ab(i) = pi*(var.id(i)^2/4);
    end
    if geom.type == "end-burner" && flag.dome && var.length(i) > geom.IL
        var.ab(i) = pi*((1 - var.height_dome(i)^2/(geom.hD^2))*(geom.OD/2)^2);
    end
    if geom.type == "rod & tube"
        if var.D1(i) <= 7.5e-2
            var.ab(i) = pi*geom.IL*(var.D1(i)+var.D2(i)) - 0.000584;
        else
            var.ab(i) = pi*geom.IL*(var.D1(i)+var.D2(i));
        end
    end 
    if geom.type == "rod & tube"
        var.ab(i) = pi*geom.IL*(var.D1(i)+var.D2(i));
    end 
    if geom.type == "bates"
        var.ab(i) = 2*pi*(var.id(i)/2)*geom.IL + 2*pi*(var.id(i)^2/4);
    end  
    % if erosive burning ---> calculate new at and eps
    if flag.nozz_erosion
        var.dt(i) = var.dt(i-1) + 2*nozz.edot*param.dt;
        var.at(i) = pi*var.dt(i)^2/4;
        var.eps(i) = nozz.ae/var.at(i);
    end

    % update burn time and iteration counter
    var.time(i) = var.time(i-1) + param.dt;
    i = i + 1;

    if i == 1e4
        a = 2+2;
    end
end

% evaluate atmospheric conditions at burnout location
[atmo.ta_b, atmo.a_b, atmo.pa_b, atmo.rho_b] = atmosisa(param.hb);

% var.pcrit = (2/(CEA.gamma+1))^(-CEA.gamma/(CEA.gamma -1));

% shutdown
while atmo.pa_b/var.pc(i-1) < param.p_ratio_burnout
    % update variables based on previous pc, tc, vol_c, and ab
    var.rb(i) = 0; % interpolated burn rate
    var.mdot_gen(i) = 0;
    var.mp_burnt(i) = 0;
    if flag.nozz_erosion
        var.mdot_out(i) = var.pc(i-1)*var.at(i-1)/CEA.cstar;
    else
        var.mdot_out(i) = var.pc(i-1)*nozz.at/CEA.cstar;
    end
    % var.mdot_sto(i) = var.mdot_sto(i-1) - var.mdot_out(i);
    var.dv_dt(i) = 0;
    var.dpc_dt(i) = CEA.R*var.tc(i-1)/var.vol_c(i-1)*(var.mdot_gen(i) - var.mdot_out(i)) - (var.pc(i-1)/var.vol_c(i-1))*var.dv_dt(i);

    % calculate performance parameters
    func = @(place) 1/place*sqrt((2/(CEA.gamma+1)*(1+(CEA.gamma-1)/2*place^2))^((CEA.gamma+1)/(CEA.gamma-1)));
    eq1 = @(place) func(place) - nozz.eps;
    perf.me(i) = fzero(eq1,1);
    perf.pe(i) = var.pc(i-1)/(1+(CEA.gamma-1)/2*perf.me(i)^2)^(CEA.gamma/(CEA.gamma-1));

    if flag.nozz_eta
        perf.cf(i) = nozz.lambda*sqrt(2*CEA.gamma^2/(CEA.gamma-1)*(2/(CEA.gamma+1))^((CEA.gamma+1)/(CEA.gamma-1))*(1-(perf.pe(i)/var.pc(i-1))^((CEA.gamma-1)/CEA.gamma))) + (perf.pe(i)-atmo.pa)/var.pc(i-1)*nozz.eps;
    else
        perf.cf(i) = sqrt(2*CEA.gamma^2/(CEA.gamma-1)*(2/(CEA.gamma+1))^((CEA.gamma+1)/(CEA.gamma-1))*(1-(perf.pe(i)/var.pc(i-1))^((CEA.gamma-1)/CEA.gamma))) + (perf.pe(i)-atmo.pa)/var.pc(i-1)*nozz.eps;
    end

    if flag.nozz_erosion
        perf.f(i)  = perf.cf(i)*var.pc(i-1)*var.at(i-1);
    else
        perf.f(i)  = perf.cf(i)*var.pc(i-1)*nozz.at;
    end

    % calculate pressure decrease
    var.pc(i) = var.pc(i-1) + var.dpc_dt(i)*param.dt;

    % calculate conditions for next interation
    var.tc(i) = CEA.tc; % need to integrate CEA for the current chamber pressure

    % calculate volume of combustion chamber
    var.vol_c(i) = var.vol_c(i-1) + var.dv_dt(i)*param.dt;

    % if erosive burning ---> calculate new at and eps
    if flag.nozz_erosion
        var.dt(i) = var.dt(i-1) + 2*nozz.edot*param.dt;
        var.at(i) = pi*var.dt(i)^2/4;
        var.eps(i) = nozz.ae/var.at(i);
    end
    % if end-burner, continue to add values to the vector
    if geom.type == "end-burner"
        var.length(i) = 0; % m
        var.id(i) = 0; % m
    end

    if flag.insulation
        var.texp(i) = 0;
        var.th_ins(i) = 0;
    end
    
    % if dome is included, continue to add values to the height of the dome
    if geom.type == "end-burner" && flag.dome 
        var.height_dome(i) = 0; % m
    end

    % update burn time and iteration counter
    var.time(i) = var.time(i-1) + param.dt;
    i = i + 1;

    if i == 135935
        a = 2+2;
    end
end


% truncate variables
var.pc = var.pc(1:(i-1),1);
var.tc = var.tc(1:(i-1),1);
var.vol_c = var.vol_c(1:(i-1),1);
if geom.type == "center-perf"
    var.id = var.id(1:(i-1),1);
    var.ab = var.ab(1:(i-1),1);
end
if flag.nozz_erosion
    var.dt = var.dt(1:(i-1),1);
    var.at = var.at(1:(i-1),1);
end
if geom.type == "end-burner"
    var.id = var.id(1:(i-1),1);
    var.ab = var.ab(1:(i-1),1);
    var.length = var.length(1:(i-1),1);
end
if geom.type == "rod & tube"
    var.D1 = var.D1(1:(i-1),1);
    var.D2 = var.D2(1:(i-1),1);
    var.ab = var.ab(1:(i-1),1);
end
if flag.insulation
    var.texp = var.texp(1:(i-1),1);
    var.th_ins = var.th_ins(1:(i-1),1);
end
if geom.type == "end-burner" && flag.dome
    var.height_dome = var.height_dome(1:(i-1),1);
end
var.rb = var.rb(1:(i-1),1);
var.mdot_gen = var.mdot_gen(1:(i-1),1);
var.mdot_out = var.mdot_out(1:(i-1),1);
% var.mdot_sto = var.mdot_sto(1:(i-1),1);
var.dv_dt = var.dv_dt(1:(i-1),1);
var.dpc_dt = var.dpc_dt(1:(i-1),1);
var.mp_burnt = var.mp_burnt(1:(i-1),1);
var.time = var.time(1:(i-1),1);
% truncate performance parameters
perf.me = perf.me(1:(i-1),1);
perf.pe = perf.pe(1:(i-1),1);
perf.cf = perf.cf(1:(i-1),1);
perf.f  = perf.f(1:(i-1),1);
end
