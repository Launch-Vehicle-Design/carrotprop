%% Solid Motor Analysis and Optimization
clear; clc; close all;

%Figure formatting preferences
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultLineMarkerSize',14);
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultFigureColor',[1,1,1]);
set(0,'DefaultTextFontSize',14);
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultTextFontName','Times-Roman');
set(0,'DefaultAxesFontName','Times-Roman');

% Volume and Mass Calculations

% griffin calcs compare
% rho_griff = 1763.66;
% mp_griff = 1465.9395;
% vol_griff = (mp_griff*(1 + residual))/rho_griff
% vol_cyl_griff = vol_griff - vol_dome
% L0_griff = vol_cyl_griff/(pi*(D0/2)^2) % m

% required initial thrust
m0 = 1780.1356;
g0 = 9.81;
t_w = 1.6;
thrust = m0*g0*t_w*0.224809

% initial parameters
prop.yi  = 1e-2*[12.00, 69.0316, 18.9684]; % propellant composition [HTPB, AP, AL]
prop.rho_i = [910 ,1950, 2710]; % propellant densities [HTPB, AP, AL]
rho = 1./sum(prop.yi./prop.rho_i); % mixture density
mp = 1278.8174; % kg
residual = 0.01;
vol_req = (mp*(1 + residual))/rho % m^3
D0 = 0.58; % m

% dome section
AR_dome = sqrt(2);
h_dome = (D0/2)/AR_dome;
vol_dome = 4/6*pi*(D0/2)*(D0/2)*(h_dome);

% cylinder section
vol_cyl = vol_req - vol_dome
L0 = vol_cyl/(pi*(D0/2)^2) % m

mp_with_resid = mp*(1 + residual)
mp_check = 1.249943446378726e+03 + vol_dome*rho

% Internal Ballistics
addpath('/Users/npixley/Documents/MATLAB/Spring 2024/AOE 4174/')

% parameters
param.temps = [-30,25,75] + 273.15; % K
param.sig_p = 0.001; % /K
param.h = 12192; % m
param.hb = 40e3; % m
param.p_ratio_burnout = 0.00554931; % from compressible aero calculator
param.g0 = 9.81; % m/s^2
param.Ru = 8.31446261815324e3;
param.dt = 1e-3;
% insulation parameters
param.tb = 130; % s
% param.tb = 300; % s
param.edot = 3*0.0000254; % mm/s
param.fs = 1.3;

% CEA outputs -> assume frozen at the throat
CEA.tc = 1598.6;
CEA.cstar = 1583.6;
CEA.gamma = 1.1980;
CEA.mw = 28.8514;
CEA.R = param.Ru/CEA.mw;

% geometry
geom.OD = D0; % m
geom.IL = 2.9; % m
geom.hD = h_dome; % m
geom.D1 = 0.2739; % m
geom.D2 = 0.3061; % m
geom.type = "rod & tube";

% propellant
prop.yi  = 1e-2*[12.00, 69.0316, 18.9684]; % propellant composition [HTPB, AP, AL]
prop.rho_i = [910 ,1950, 2710]; % propellant densities [HTPB, AP, AL]
prop.rho = 1./sum(prop.yi./prop.rho_i); % mixture density
% factor = 2.45; % multiplication factor to increase the burn rate
% prop.rb = 1e-3*[7.84, 8.82, 10.05, 10.62, 11.29]*factor; % m/s -> burn rate [HTPB, AP, AL]
% prop.pc = [3, 5, 7, 9, 11]; % MPa -> chamber pressure [HTPB, AP, AL]
prop.rb = 0.107*[0.00721, 0.0072, 0.00722, ...
                     0.01372, 0.01373, 0.01378, ...
                     0.01812, 0.018, 0.0182, ...
                     0.0213, 0.023, 0.0214, ...
                     0.0239, 0.0239, 0.0241, ...
                     0.0291, 0.0286, 0.0276];
prop.pc = 1.2*[0.689, 0.687, 0.692, ...
                 3.446, 3.444, 3.45, ...
                 6.893, 6.873, 6.893, ...
                 10.339, 10.342, 10.34, ...
                 13.785, 13.79, 13.787, ...
                 20.684, 20.7, 20.65];

% solve for a and n
prop.pc_log = log(prop.pc);
prop.rb_log = log(prop.rb);
p = polyfit(prop.pc_log,prop.rb_log,1);
prop.a = exp(p(2));
prop.n = p(1);

% plot burning rate data
figure()
hold on
grid on
scatter(prop.pc,prop.rb,150,"black",'LineWidth',1.5)
set(gca,'xscale','log')
set(gca,'yscale','log')
x = linspace(prop.pc(1),prop.pc(end));
y = prop.a*x.^prop.n;
plot(x,y)
xlabel('$P_c$ [MPa]')
ylabel('$r_b$ [m/s]')
title('Burning Rate Regression Data - $r_b = a*{P_c}^n$')
hold off

% atmospheric conditions
[atmo.ta, atmo.a, atmo.pa, atmo.rho] = atmosisa(param.h);

% define nozzle characteristics
nozz.dt = 0.1; % m -> diameter of the throat
nozz.at = pi*nozz.dt^2/4; % Throat Area (m^2)
nozz.eps = 18.75; % expansion ratio
nozz.ae = nozz.eps*nozz.at; % Exit Area (m^2)
nozz.edot = 5e-5; % m/s
nozz.alpha = 15; % degrees
nozz.lambda = 1/2*(1 + cosd(nozz.alpha)); % correction factor

% define ignition parameters
ignit.pc(1) = 101325*2.2; % Pa  -> pressure after ignition
ignit.tc(1) = CEA.tc; % K   -> temperature after ignition
ignit.vol_c(1) = .0361; % m^3 -> initial volume of motor casing

% [ambient -> no nozzle erosion]
flag = struct('nozz_erosion', false, ...
              'temp_sens_high', false, ...
              'temp_sens_low', false, ...
              'temp_sens_ambient', true, ...
              'nozz_eta', true,...
              'dome', true,...
              'insulation', false);

% calculate performance
tic
[prop,nozz,var,perf] = internalballistics_V1(param,CEA,prop,geom,nozz,ignit,flag)
toc

% % [ambient -> include nozzle erosion]
% flag = struct('nozz_erosion', true, ...
%               'temp_sens_high', false, ...
%               'temp_sens_low', false, ...
%               'temp_sens_ambient', true, ...
%               'nozz_eta', false,...
%               'dome', false);


% % calculate performance
% tic
% [prop_e,nozz_e,var_e,perf_e] = internalballistics(param,CEA,prop,geom,nozz,ignit,flag)
% toc

% chamber pressure vs. time
figure()
hold on
grid on
plot(var.time(:,1),var.pc(:,1)*1e-6,'LineWidth',3)
% plot(var_e.time(:,1),var_e.pc(:,1)*1e-6,'LineWidth',3)
xlabel('Time [s]','FontSize',25)
ylabel('$P_c$ [MPa]','FontSize',25)
title('Chamber Pressure','FontSize',25)
% legend('Non-Erosive','Erosive','Interpreter','latex','Location','northwest','FontSize',15)
hold off

% thrust vs. time
figure()
hold on
grid on
plot(var.time(:,1),perf.f(:,1)*0.224809,'LineWidth',3)
% plot(var_e.time(:,1),perf_e.f(:,1)*0.224809,'LineWidth',3)
xlabel('Time [s]','FontSize',40)
ylabel('Thrust [lb$_f$]','FontSize',40)
title('Thrust','FontSize',40)
% legend('Non-Erosive','Erosive','Interpreter','latex','Location','northwest','FontSize',15)
hold off

% mass flow vs. time
figure()
hold on
grid on
plot(var.time(:,1),var.mdot_out(:,1)*2.20462,'LineWidth',3)
% plot(var_e.time(:,1),var_e.mdot_out*2.20462,'LineWidth',3)
xlabel('Time [s]','FontSize',40)
ylabel('$\dot{m}$ [$\frac{lb_m}{s}$]','FontSize',40)
title('Mass Flow','FontSize',40)
% legend('Non-Erosive','Erosive','Interpreter','latex','Location','northwest','FontSize',15)
hold off

% mass flow vs. time
figure()
hold on
grid on
plot(var.time(1:(end-259),1),var.mdot_gen(1:(end-259),1)*2.205,'LineWidth',3)
plot(var.time(:,1),var.mdot_out(:,1)*2.205,'LineWidth',3)
% plot(var.time(:,1),var.mdot_sto(:,1),'LineWidth',1.5)
xlabel('Time [s]','FontSize',40)
ylabel('$\dot{m}$ [$\frac{lb_m}{s}$]','FontSize',40)
title('Mass Flow','FontSize',40)
legend('$\dot{m}_{gen}$','$\dot{m}_{out}$','Interpreter','latex','Location','southwest','FontSize',40)
hold off

% display relevant outputs
% % % need to calculate this at the end to verify a residence time of ~ 1 ms
% delta_t_estimate = 1/5*(7.5e6*0.7317)/((1-0.614)*(Ru/24.4184*2000*15));

% % calculate insulation mass
% rho_insulation = 76.1; % lb/ft^3
% rho_insulation = rho_insulation*16.018; % kg/m^3
% Y = var.length;
% X = 2*pi*var.th_ins;
% vol_ins = trapz(Y,X) % m^3
% 
% disp('Insulation Mass [lbm]')
% insul_mass = vol_ins*rho_insulation*2.20462; % lbm
% disp(insul_mass)

% define propellant characteristics/constraints
% prop.mp = 1460.9017; % kg [ original propellant mass]
% prop.mp = 1437.2538; % kg

% burn rate (avg)
disp('rb_avg [in/s]')
rb_avg = mean(var.rb)*1000*0.0393701;
disp(rb_avg)
disp('rb_avg [mm/s]')
rb_avg_metric = mean(var.rb)*1e3;
disp(rb_avg_metric)

% contraction ratio
disp('eps_c')
eps_c = pi*(geom.OD/2)^2/(nozz.at);
disp(eps_c)

% mass flow
disp('mdot_avg [lbm/s]')
mdot = (mean(var.pc))*nozz.at/CEA.cstar*2.20462;
disp(mdot)

% ignitor propellant mass calculations
disp('mp_ignitor [kg]')
mp_ign = 0.12*(pi*(geom.D2^2/4 - geom.D1^2/4)*geom.IL)^0.7;
disp(mp_ign)

% insulation thickness calculations
t_b = max(var.time(:,1))
edot_case = 3; % mil/s
edot_case = edot_case*0.0254; % mm
f_s = 1.3;
d_max = t_b*edot_case*f_s; % mm

disp('Max Thickness [in]')
d_max_in = d_max*0.0393701; % in
disp(d_max_in)

% area of rod support
width = 0.075;
thickness = 0.004; 
area = 2*(width*thickness) - thickness*thickness;
disp('Support Cross-Sectional Area [m^2]')
disp(area)

% 7.5 cm * 0.4 cm
% 7.5 cm * 0.4 cm
% - 0.4 cm * 0.4 cm
% 
% 7.5*0.4*2 - 0.4*0.4=5.84 cm^2
% 
% L0 = 2.9 m - 290 cm

% V0 = 290*5.84=1,693.6 cm^3


%% % % will need to include starting and ending atmospheric conditions
% linerally space the data points between these values to approximate

% thrust coefficient

% thrust vs. time
figure()
hold on
grid on
plot(var.time(1:(end-64),1),perf.f(1:(end-64),1)*0.224809,'LineWidth',3)
xlabel('Time [s]','FontSize',25)
ylabel('$\tau [lb_f]$','FontSize',25)
title('Thrust Profile','FontSize',25)
hold off

% mass flow vs. time
figure()
hold on
grid on
plot(var.time(1:(end-64),1),var.mdot_gen(1:(end-64),1)*2.205,'LineWidth',3)
plot(var.time(:,1),var.mdot_out(:,1)*2.205,'LineWidth',3)
% plot(var.time(:,1),var.mdot_sto(:,1),'LineWidth',1.5)
xlabel('Time [s]','FontSize',25)
ylabel('$\dot{m}$ [$\frac{lb_m}{s}$]','FontSize',25)
title('Mass Flow','FontSize',25)
legend('$\dot{m}_{gen}$','$\dot{m}_{out}$','$\dot{m}_{stored}$','Interpreter','latex','Location','southwest','FontSize',25)
hold off

% chamber pressure vs. time
figure()
hold on
grid on
plot(var.time(:,1),var.pc(:,1)*9.8692e-6,'LineWidth',3)
xlabel('Time [s]','FontSize',25)
ylabel('$P_c$ [atm]','FontSize',25)
title('Chamber Pressure','FontSize',25)
hold off
toc

% %% Animated 2D Gif of solid propellant burn
% % Define the number of frames
% numFrames = 100;
% 
% % Pull length of fuel grain
% [M,I] = max(var.length);
% 
% % Initial size and position of the rectangle
% width = max(M);
% initialHeight = prop.d0;
% positionX = 0.0; % X-coordinate of the lower left corner
% positionY = 0.0; % Y-coordinate of the lower left corner
% 
% % Initialize the figure
% figure;
% 
% % Loop to generate each frame
% for k = 1:I
%     % Calculate the new size of the rectangle
%     width = width - 1000*var.rb(k)*dt;
%     if width <= 0
%         width = 0;
%     end
%     height = prop.d0;
% 
%     % Clear the figure to update the rectangle in each iteration
%     clf;
% 
%     % Plot the shrinking rectangle
%     rectangle('Position', [positionX positionY width height], 'EdgeColor', 'r', 'LineWidth', 2);
%     axis([0 3 0 2]); % Fix the axis so the shrinking effect is visible
%     axis equal; % Keep the aspect ratio of the plot consistent
% 
%     % Capture the frame
%     frame = getframe(gcf);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
% 
%     % Write to the GIF File 
%     if k == 1
%         imwrite(imind,cm,'shrinkingRectangle.gif','gif', 'Loopcount',inf);
%     else
%         imwrite(imind,cm,'shrinkingRectangle.gif','gif','WriteMode','append');
%     end
% 
%     if width == 0
%         break
%     end
% end

% %% pre fire
% close all;
% 
% % Load the images
% img1 = imread('solidmotor_after.jpg'); % Background image
% img2 = imread('solidmotor_before.jpg'); % Image to reveal
% 
% % Ensure both images are the same size
% assert(all(size(img1) == size(img2)), 'Images must be the same size.');
% 
% % Number of frames for the animation
% numFrames = 150;
% 
% % Initial dimensions for the shrinking rectangle
% rectWidth = size(img1, 2); % Full width
% rectHeight = size(img1, 1); % Full height
% centerX = size(img1, 2) / 2;
% centerY = size(img1, 1) / 2;
% 
% % Initialize the figure
% figure;
% 
% % Animation loop
% for k = 1:numFrames
%     % Calculate new dimensions for the shrinking rectangle
%     % Here, only xStart changes to simulate the shrinking from the left side.
%     width = max(1, rectWidth * (1 - k / numFrames));
% 
%     % xStart moves to the right as the rectangle shrinks
%     xStart = max(1, round(rectWidth - width));
%     % xEnd stays fixed at the full width of the image
%     xEnd = size(img1, 2);
% 
%     height = rectHeight;  % The height remains constant
% 
%     % Create a binary mask for the current frame
%     mask = false(size(img1, 1), size(img1, 2));
%     yStart = max(1, round(centerY - height / 2));
%     yEnd = min(size(img1, 1), round(centerY + height / 2));
%     mask(yStart:yEnd, xStart:xEnd) = true;
% 
%     % Blend the images based on the mask
%     compositeImg = img1;
%     for c = 1:size(img1, 3) % For each color channel
%         channel = compositeImg(:,:,c);
%         channel(mask) = img2(yStart:yEnd, xStart:xEnd, c);
%         compositeImg(:,:,c) = channel;
%     end
% 
%     % Display the composite image
%     imshow(compositeImg);
% 
%     % Capture the frame
%     frame = getframe(gcf);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
% 
%     % Write to the GIF file
%     if k == 1
%         imwrite(imind, cm, 'revealingImage.gif', 'gif' ,'DelayTime',0.000001, 'Loopcount', inf);
%     else
%         imwrite(imind, cm, 'revealingImage.gif', 'gif' ,'DelayTime',0.000001, 'WriteMode', 'append');
%     end
% end
% 
% %% post fire
% close all;
% 
% % Load the images
% img1 = imread('solidmotor_after.jpg'); % Background image
% img2 = imread('solidmotor_before.jpg'); % Image to reveal
% 
% % Load the flame image with its alpha channel for transparency
% [flameImg] = imread('idk_test_123.png');
% % flameImg = imrotate( flameImg , 90 );
% 
% % Ensure all images are the same height
% assert(size(img1, 1) == size(img2, 1), 'The first and second images must be the same height.');
% assert(size(img1, 1) == size(flameImg, 1), 'The flame image must match the height of the other images.');
% 
% % Number of frames for the animation
% numFrames = 150;
% 
% % Initial dimensions for the shrinking rectangle
% rectWidth = size(img1, 2); % Full width
% rectHeight = size(img1, 1); % Full height
% 
% % Initialize the figure
% figure;
% 
% % Flame image properties
% flameWidth = size(flameImg, 2);
% flameHeight = size(flameImg, 1);
% 
% % Animation loop
% for k = 1:numFrames
%     % Calculate new dimensions for the shrinking rectangle
%     width = max(1, rectWidth * (1 - k / numFrames));
% 
%     % xStart moves to the right as the rectangle shrinks
%     xStart = max(1, round(rectWidth - width));
%     xEnd = size(img1, 2);
% 
%     % Create the mask for the current frame
%     mask = true(size(img1, 1), size(img1, 2));
%     mask(:, xStart:xEnd) = false;
% 
%     % Blend the images based on the mask
%     compositeImg = img1;
%     compositeImg(repmat(~mask, [1 1 3])) = img2(repmat(~mask, [1 1 3]));
% 
%     % Calculate the position for the flame image
%     flamePosX = max(1, xStart - flameWidth);
% 
%     % Place the flame image on top of the background and the receding rectangle
%     % Note: This assumes the flame image background is black and we're using that
%     % as a way to not show the flame where it's black.
%     for c = 1:3
%         compositeChannel = compositeImg(1:flameHeight, flamePosX:flamePosX+flameWidth-1, c);
%         flameChannel = flameImg(:, :, c);
%         blackPixels = flameChannel == 0; % Find the black background pixels
%         flameChannel(blackPixels) = compositeChannel(blackPixels); % Keep original background where the flame is black
%         compositeImg(1:flameHeight, flamePosX:flamePosX+flameWidth-1, c) = flameChannel;
%     end
% 
%     % Display the composite image
%     imshow(compositeImg);
% 
%     % Capture the frame
%     frame = getframe(gcf);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
% 
%     % Write to the GIF file
%     if k == 1
%         imwrite(imind, cm, 'revealingImageWithFlame.gif', 'gif', 'Loopcount', inf);
%     else
%         imwrite(imind, cm, 'revealingImageWithFlame.gif', 'gif', 'WriteMode', 'append');
%     end
% end


%% functions
% burning area functions
function area = area_endburn(d0)
    area = pi*d0^2/4;
end
