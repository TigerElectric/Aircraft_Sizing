%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bernardo Pacini & Leif Fredericks                                     %%
%% Hybrid-Electric General Aviation Aircraft                             %%
%% Preliminary Design Calculations                                       %%
%% Created: 1/23/2018                                                    %%
%% Modified: 2/15/2018                                                   %%
%%                                                                       %%
%% Description: This code will output preliminary aircraft design        %%
%% calculations to a .txt file and a carpet plot with constraints        %%
%%                                                                       %%
%% Dependencies:| aircraft_weight.m | aircraft_carpetplot.m | Atmos.m |  %%
%%              | calculate_alpha.m | calculate_beta.m | TS_converter.m |%% 
%%              | convert_to_imperial.m | dynamic_pressure.m |           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Trial_Name  = 'Trial 1' ;
Description = ''        ;
% ON/OFF SELECTORS
output      = 0         ; %1/0 for output files
carpet_plot = 1         ; %1/0 for carpet plot
weight      = 1         ; %1/0 for weight
graph       = 0         ; %1/0 for graph of weight equation solution

% PERFORMANCE ESTIMATES
AR          = 7         ; % Aspect ratio: ESTIMATE then ITERATE
LDmax       = AR+10     ;  % early guess, should change **********************
LDC         = 15        ; % L/D at cruise conditions (based on design C_L),
LD          = 18        ; % L/D max: ESTIMATE (LDC / 0.94) then ITERATE
eta         = 0.8       ; % propulsion efficiency
e = 1.78*(1-0.045*AR^(0.68)) - 0.64;
k = 1/(pi*AR*e); % induced drag correction factor, 

ws          = 75        ; % wing loading (to check design point Landing FL)

% AIRCRAFT REQUIREMENTS 
altitude_ci = 20000     ; %cruise altitude, ft
altitude_fi = 0000      ; %airfield alitude, ft

% FLIGHT REQUIREMENTS
V_cruise    = 200       ;
V_stall     = 61        ; % knots FAR 23
V_approach  = V_stall*1.3; %knots
L_takeoff   = 1800      ; %ft REQUIREMENT
L_landing   = 1800      ; %ft REQUIREMENT
V_climb     = V_stall*1.2; %for now (see aircraft_mass.m)
rate_climb  = 1300      ; %ft/min
altitude_climbi = 0     ; %ft, for now (see aircraft_mass.m)
theta_app   = 3.04      ; %approach angle, deg
n_max       = 6         ; % maximum load factor (FAR 23)
n_min       = -0.5      ; % minimum load factor (FAR 23)

R           =   750;    % [nmi]     RFP ***********************************
t_emergency_guess =   60 / 60;% [Hrs]     FAR-ish *******************************
passengers  =   5;      %           RFP ***********************************
crew        =   1;      %           RFP ***********************************


% CRUISE PARAMETERS
gamma       = 1.4       ; % specific heat ratio cp/cv, for air
g           = 32.174    ; %ft/s^2

Clmax_to    = 1.50      ; % assumed: ESTIMATE then ITERATE
WP_guess    = 15        ; 
WS_guess    = 26        ;
C_D0        = 0.03      ; % assumed (at cruise): ESTIMATE then ITERATE
TW          = 550*eta/V_climb / WP_guess; % to T/W, can use WP cruise cause what engine sees
Clmax_climb = Clmax_to/1.44; % Hamburg ch. 5
Clmax_land  = 2.10      ; % assumed: ESTIMATE then ITERATE

% CARPET PLOT LIMITS
carpet_x_lim = [1 100]  ; % W/S ratio, generally between 50 and 100
carpet_y_lim = [0 50]   ; % T/W ratio, generally between .4 and .7

SHOW      =   0;      % SWITCH to show design points on non-log Ragone
%%  STARRED ***** should be fed in as function from Prelim_Sizing
% Performance Parameters
MTOW        =   12000   ;   % FAR is 12,500 lb for FAR 23 classification*******
bsfc        =   0.3     ;    % Aproximate based on Leif's engine # *************

weight_max  =   10^6;   % approximate maximum weight for iterations *******
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%DO NOT MODIFY BELOW THIS POINT%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTION CALLING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Surface Area
if carpet_plot == 1
%try
[HPW_takeoff_ws, HPW_climb_ws, HPW_cruise_ws, HPW_turning_max_ws, ...
    HPW_turning_min_ws, HPW_endurance_ws, W_S] = ...
    aircraft_carpetplot_power(V_cruise, AR, eta, WP_guess, WS_guess, ...
    altitude_ci, altitude_fi,...
    altitude_climbi, V_approach, V_stall, Clmax_to, Clmax_climb, ...
    Clmax_land, L_takeoff, L_landing, V_climb, rate_climb, theta_app,...
    C_D0, gamma, g, carpet_x_lim, ...
    carpet_y_lim, LD, LDC, n_max, n_min, k, e);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WEIGHT CALCULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if weight == 1
% N_prop for loiter is hardcoded to 0.7 as per Raymer, but may need to change
% 190 lb passengers and 30 lb for luggage

% MTOW[lb], R[nmi], V_cruise[Kts], V_stall[Kts], L_takeoff[ft],
% fpm_climb[ft/min], cruise_alt[ft], t_emergency_guess[Hrs]


[W_TO, W_fuel, W_batt, W_empty, W_pay, W_batt_TO, W_batt_CM, W_batt_EM, t_EM]...
    = Weight_est(MTOW, bsfc,...
AR, C_D0, Clmax_to, LDmax, eta, HPW_cruise_ws, HPW_climb_ws, ...
HPW_takeoff_ws, HPW_turning_min_ws, W_S,...
weight_max, R, V_cruise, V_stall, L_takeoff, rate_climb, altitude_ci,...
t_emergency_guess, passengers, crew, SHOW)

t_EM/60*V_cruise
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TAKEOFF DISTANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% try
% if TO_landing_lengths == 1
% [ TOFL , index] = TO_distance(V_stall, Clmax_to, s_ref, ... 
%     Thrust, W_TO, Number_Engines, altitude_fi);
% disp(sprintf('%0.0f Required Takeoff Eield Length (ft)', TOFL));
% end
% catch
%     errordlg('Please enable both weight and Carpet Plot calculations')
%     error('Please enable both weight and Carpet Plot calculations')
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT TO TEXT FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if output == 1
folder_name = uigetdir('C:\','Select Working Directory');
mkdir(sprintf('%s/%s', folder_name, Trial_Name));

text = sprintf('Preliminary Aircraft Design Calculations - %s.txt',...
    Trial_Name);
fid = fopen(sprintf('%s/%s/%s',folder_name, Trial_Name, text),'w');
fprintf(fid, sprintf('Preliminary Aircraft Design Calculations- %s \n',...
    Trial_Name));
fprintf(fid, sprintf('Description: %s \n \n',Description));
if weight == 1
fprintf(fid, sprintf('%0.0f Takeoff Weight (lbm) \n', W_TO)); 
fprintf(fid, sprintf('%0.0f Fuel Weight (lbm) \n', W_fuel));
fprintf(fid, sprintf('%0.0f Empty Weight (lbm) \n \n', W_empty));
end
if carpet_plot == 1
fprintf(fid, sprintf('%0.0f Thrust (lbf) \n', thrust));
fprintf(fid, sprintf('%0.0f Reference Area (ft^2) \n \n', s_ref));
end
if TO_landing_lengths == 1
fprintf(fid, sprintf('%0.0f Takeoff Field Length (ft) \n \n',TOFL));
fprintf(fid,sprintf('%0.0f Landing Field Length (ft) \n \n',Landing_Dist));
end

fprintf(fid, sprintf('------------------------------------------------'));
fprintf(fid, sprintf('\nInput Parameters: \n \n'));

fprintf(fid, sprintf('%0.2f Cruise Mach Number\n', V_cruise)); 
fprintf(fid, sprintf('%0.0f Range (Nm)\n', R)); 
fprintf(fid, sprintf('%0.1f Aspect Ratio\n', AR)); 
fprintf(fid, sprintf('%0.2f Oswald Efficiency Factor\n', e)); 
fprintf(fid, sprintf('%0.3f TSFC\n', tsfc)); 
fprintf(fid, sprintf('%0.0f Altitude_ci (ft) \n', altitude_ci)); 
fprintf(fid, sprintf('%0.0f Altitude_fi (ft) \n', altitude_fi)); 
fprintf(fid, sprintf('%0.0f Passengers \n', passengers)); 
fprintf(fid, sprintf('%0.0f Crew \n', crew)); 
fprintf(fid, sprintf('%0.0f Baggage (lbm) \n', baggage)); 
fprintf(fid, sprintf('%0.2f Loiter Duration (hrs) \n \n', loiter_dur)); 

fprintf(fid, sprintf('%0.1f Approach Velocity (knots) \n', V_approach)); 
fprintf(fid, sprintf('%0.2f Stall Velocity (knots) \n', V_stall)); 
fprintf(fid, sprintf('%0.2f Takeoff CL_max \n', Clmax_to)); 
fprintf(fid, sprintf('%0.2f Land CL_max \n', Clmax_land)); 
fprintf(fid, sprintf('%0.2f Climb Mach Number \n', V_climb)); 
fprintf(fid, sprintf('%0.2f Rate Climb (ft/min) \n', rate_climb)); 
fprintf(fid, sprintf('%0.2f Climb Altitude (ft) \n', altitude_climbi)); 
fprintf(fid, sprintf('%0.2f Approach angle [theta] (deg)\n\n', theta_app)); 

fprintf(fid, sprintf('%0.3f Coefficient of Drag \n', C_D0)); 
fprintf(fid, sprintf('%0.3f Coefficient of Residual Drag \n', C_DR_c)); 
fprintf(fid, sprintf('%0.3f Induced Drag Correction [K1] \n', K1_c)); 
fprintf(fid, sprintf('%0.3f Viscous Drag Correction [K2] \n', K2_c)); 
fprintf(fid, sprintf('%0.3f Specific Heat Ratio [gamma] \n', gamma)); 
fprintf(fid, sprintf('%0.1f TR \n', TR)); 
fprintf(fid, sprintf('%0.1f Gravitational Constant (ft/s^2) \n', g)); 

if carpet_plot == 1
print(sprintf('%s/%s/Carpet Plot %s', folder_name, Trial_Name, ...
Trial_Name),'-dpng')
end
end