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
addpath('/Users/bernardo_pacini/Documents/MAE 340D-HEGAA')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prompt = {'File Name'};
dlg_title = 'User Input';
num_lines = 1;
defaultans = {'Filename'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
Trial_Name  = answer{1} ;
Description = ''        ;
% ON/OFF SELECTORS
output      = 1         ; %1/0 for output files
carpet_plot = 1         ; %1/0 for carpet plot
weight      = 1         ; %1/0 for weight
SHOW        = 0         ; %1/0 for graph of weight equation solution and Ragone
%%
% PERFORMANCE ESTIMATES (CHANGE EVERY TIME!!!!!!!!!!!!!!!)
LDC         = 12.9        ; % L/D at cruise conditions (based on AoA and CL)
C_D0        = 0.020      ; % assumed (at cruise): ESTIMATE then ITERAT
LDL         = 15.4        ; % L/D at loiter conditions (based on AoA and CL)
Cl_stall    = 1.3      ; % Probably at 15 deg, need to find a source, or use airfoil NACA 63-415
                          % calculate at landing(stalling) Mach Number
bsfc        = 0.30      ; % Aproximate based on Leif's engine # 
 
% SECOND LEVEL DESIGN PARAMETERS (May Change Slightly)
AR          = 7         ; % Aspect ratio
eta         = 0.8       ; % Propeller efficiency
TOP = 200; % Raymer fig 5.4 for 1800 ish with 50 obstacle
emergency_range = 105; % Set emergency range in nmi
PDENS_HP    =   [3000; 440; 300]*0.000608277;
EDENS_HP    =   [120; 220; 300]*0.000608277;


% DERIVED PARAMETERS (Shouldn't Change Much)
e = 1.78*(1-0.045*AR^(0.68)) - 0.64;
k = 1/(pi*AR*e); % induced drag correction factor, 

Clmax_land  = Cl_stall+1.35; % http://www.fzt.haw-hamburg.de/pers/Scholz/
                             % HOOU/AircraftDesign_5_PreliminarySizing.pdf
                             % Need to use a fowler or double-slotted
Clmax_to    = 0.8*Clmax_land;   % Last time we said this is in Raymer
Clmax_climb = Clmax_to/1.44;


% SIZING PARAMETERS (Probably won't change much)
AR_v            =   1.6; % Aspect Ratio of Vertical Tail, Raymer 113 4.5.4
AR_h            =   5;   % Aspect Ratio of Horizontal Tail  
lambda          =   0.6; % Taper Ratio of Wings,per Raymer 83 i am trying 0.6 cause it looks too weird with .4
lambda_v        =   0.4; % Taper Ratio of vertical Tail, Raymer 113
lambda_h        =   0.6; % Taper Ratio of Horizontal Tail, Raymer 113
sweep_LE        =   0;   % Sweep [deg] based on M 0.33 max. Raymer 80 4.3.2
sweep_LE_v      =   35;  % Sweep [deg] of Vertical Tail, Raymer 112
sweep_LE_h      =   5;   % Sweep [deg] of Horizontal Tail, Raymer 112
c_VT            =   0.04; % VT volume coeffieient Initial Guess: 0.09 per Raymer pg. 160 and Martinelli
c_HT            =   0.70; % HT volume coeff. Raymer 6.5.3
dihedral        =   5; % rad, guess from Raymer Table 4.2 for low-wing placement, 
upsweep         =   25;
D_fus           =   5;   % Easier to set, just make sure not too diff from recommended
a               =   4.37; % Fuselage Constant
C               =   0.23; % Fuselage Constant
fineness        =   7; % Raymer 6.5.1

% AIRCRAFT REQUIREMENTS (Shouldn't Change Much)
altitude_ci = 20000     ; %cruise altitude, ft
altitude_fi = 0000      ; %airfield alitude, ft
altitude_li = 5000      ; % Loiter Altitude, ft (can set wehrever cause turbo)
MTOW        = 12000;   % FAR is 12,500 lb for FAR 23 classification*******
weight_max  = 10^15;   % approximate maximum weight for iterations *******

% FLIGHT REQUIREMENTS (Change Between Big and Small)
V_cruise        = 200           ; % This is Mach 0.3 at 20,000
V_stall         = 61            ; % knots FAR 23
V_approach      = V_stall*1.3   ; %knots
L_takeoff       = 1800          ; %ft REQUIREMENT
L_landing       = 1800          ; %ft REQUIREMENT
V_climb         = V_stall*1.2   ; %for now (see aircraft_mass.m)
rate_climb      = 1300          ; %ft/min
loiter_dur      = 0.75          ; % [hr] Pretty sure this is in Raymer
altitude_climbi = 0             ; %ft, for now (see aircraft_mass.m)
n_max           = 4.4           ; % maximum load factor (FAR 23)
n_min           = -0.4*n_max    ; % minimum load factor (FAR 23)
R               = 750           ; % [nmi] RFP
passengers      = 5             ; % RFP
crew            = 1             ; % RFP 

% CRUISE PARAMETERS (Don't Change)
gamma       = 1.4       ; % specific heat ratio cp/cv, for air
g           = 32.174    ; %ft/s^2

% CARPET PLOT LIMITS (Shouldn't Need To Change)
carpet_x_lim = [1 100] ; % W/S ratio, generally between 50 and 100
carpet_y_lim = [0 50]    ; % T/W ratio, generally between .4 and .7

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%DO NOT MODIFY BELOW THIS POINT%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTION CALLING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Surface Area
[HPW_takeoff_ws, HPW_climb_ws, HPW_cruise_ws, HPW_turning_max_ws, ...
    HPW_turning_min_ws, HPW_endurance_ws, W_S, HPW_match] = ...
    aircraft_carpetplot_power(V_cruise, AR, eta, ...
    altitude_ci, altitude_fi,...
    altitude_climbi,altitude_li, V_approach, V_stall, Clmax_to, Clmax_climb, ...
    Clmax_land, L_takeoff, L_landing, V_climb, rate_climb,...
    C_D0, g, carpet_x_lim, carpet_y_lim, n_max, n_min,LDL,TOP);


[airDens_ci, airPres_ci, temp_ci, soundSpeed_ci] = Atmos(altitude_ci*0.3048);
[airDens_li, airPres_li, temp_li, soundSpeed_li] = Atmos(altitude_li*0.3048);
rho_air_ci       =              airDens_ci*0.062428;    % lb/ft^3
rho_air_li       =              airDens_li*0.062428;    % lb/ft^3
c_air_ci         =              soundSpeed_ci*1.94384;  % kts
c_air_li         =              soundSpeed_li*1.94384;  % kts
V_cruise_fps     =              V_cruise*1.68781;       % ft/s
[airDens_fi, airPres_fi, temp_fi, soundSpeed_fi] = Atmos(altitude_fi*0.3048);
rho_air_fi       =              airDens_fi*0.062428;    % lb/ft^3
c_air_fi         =              soundSpeed_fi*1.94384;  % kts
V_TO_fps         =              1.2*V_stall*1.68781;    % ft/s

CL_cruise       = W_S * g / (0.5 * rho_air_ci * V_cruise_fps^2);
M_cruise        = V_cruise/c_air_ci;
M_land          = V_stall/c_air_fi;
V_loiter_fps    = sqrt(2*W_S*g/rho_air_ci * sqrt(k / (3*C_D0 )));
V_loiter        = sqrt(2*W_S*g/rho_air_ci * sqrt(k / (3*C_D0 )))/1.68781; %kts
M_loiter        = V_loiter/c_air_li;
CL_loiter       = W_S * g / (0.5 * rho_air_li * V_loiter_fps^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WEIGHT CALCULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% N_prop for loiter is hardcoded to 0.7 as per Raymer, but may need to change
% 190 lb passengers and 30 lb for luggage

% MTOW[lb], R[nmi], V_cruise[Kts], V_stall[Kts], L_takeoff[ft],
% fpm_climb[ft/min], cruise_alt[ft], t_emergency_guess[Hrs]
[W_TO, W_fuel, W_batt, W_empty, W_pay, W_batt_TO, W_batt_CM, W_batt_EM, t_EM, TBat]...
    = Weight_est_fixbatt(MTOW, bsfc,...
AR, C_D0, Clmax_to, LDC, LDL, eta, HPW_cruise_ws, HPW_climb_ws, HPW_takeoff_ws, HPW_turning_min_ws, W_S,...
weight_max, R, V_cruise, V_stall, L_takeoff, rate_climb, altitude_ci,...
emergency_range, passengers, crew,loiter_dur, SHOW, PDENS_HP, EDENS_HP, V_loiter,HPW_endurance_ws);

Parameters_W = {'MTOW';'Fuel Weight';'Total Battery Weight';'Empty Weight';...
    'Payload Weight';'TO Battery Weight';'Climb Battery Weight';...
    'Emergency Battery Weight'};
Weights_lb = [W_TO; W_fuel; W_batt; W_empty; W_pay; W_batt_TO; W_batt_CM;...
    W_batt_EM];
W = table(Weights_lb, 'RowNames', Parameters_W);
disp(W);

Emergency_Dist  = t_EM/60*V_loiter;
S               = W_TO / W_S;
Power           = W_TO / HPW_match;


Parameters_Y = {'Match W/S'; 'Match W/P'; 'Surface Area (ft)';...
    'Power Requirement (HP)'; 'Mach Landing'; 'Mach Cruise';...
    'Mach Loiter';'CL at Cruise'; 'CL at Loiter'; 'Loiter Velocity';...
    'Emergency Time (min)'; 'Emergency Distance (nmi)'};
Defining_Properties = [W_S; HPW_match; S; Power; M_land; M_cruise; ...
    M_loiter;CL_cruise; CL_loiter; V_loiter; t_EM; Emergency_Dist];
Y = table(Defining_Properties, 'RowNames', Parameters_Y);
disp(Y);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GEOMETRIC SIZING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORTANT SIZING NOTES ON THIS
% The fuselage rear end point is assumed to be 1/4 of the diameter above
% the center line, i.e. 75% of the diameter above the point where it must
% start to angle up. Set the last point in the fuselage to +(1/4)*D_fus in
% the z axis in VSP
% These are fuselage parameters in ft for a single engine GA 6.5.1 Raymer
L_engine = 0.0083325 * Power; % [ft] See thoughts and sources in sizing folder
L_nose = 2*L_engine; % [ft] Approximate room for engine mounting, no prop included
% FUSELAGE SIZING, no prop cone included
[L_fus, L_HT, L_VT, V_fus, D_fus_rec, S_wet_fus, L_fus_angled,V_batt] = ...
    hybrid_fuselage(a, W_TO, W_batt, C, fineness, L_nose, upsweep, D_fus);
R_fus = D_fus/2; % Radius of fuselage for wing surface placements
Fuselage_Parameter = [L_fus;D_fus;D_fus_rec;V_fus;S_wet_fus;L_fus_angled;...
                                L_fus-L_fus_angled; L_nose; L_nose];
fus_param_names = {'Length','Max Diameter','Recommended Diameter'...
    ,'Volume Fuselage','Wet Area','Angled Length', 'Angled Starts Maximum',...
    'Nose Length','Wing Position (Oversized Tail)'};
F = table(Fuselage_Parameter, 'RowNames', fus_param_names);
disp(F);

[TSize] = hybrid_sizing(S, AR, lambda, sweep_LE, AR_v, AR_h, lambda_v, ...
    lambda_h, sweep_LE_v, sweep_LE_h, c_VT, c_HT, L_VT, dihedral,...
    L_fus, D_fus);
disp(TSize)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

fprintf(fid, sprintf('%0.0f Takeoff Weight (lbm) \n', W_TO)); 
fprintf(fid, sprintf('%0.0f Fuel Weight (lbm) \n', W_fuel));
fprintf(fid, sprintf('%0.0f Battery Weight (lbm) \n', W_batt));
fprintf(fid, sprintf('%0.0f Empty Weight (lbm) \n', W_empty));
fprintf(fid, sprintf('%0.0f Payload Weight (lbm) \n', W_pay));
fprintf(fid, sprintf('%0.1f Take-Off Battery Weight (lbm) \n', W_batt_TO));
fprintf(fid, sprintf('%0.0f Climb Battery Weight (lbm) \n', W_batt_CM));
fprintf(fid, sprintf('%0.0f Emergency Battery Weight (lbm) \n\n', W_batt_EM));

fprintf(fid, sprintf('%0.2f Match Point W/S \n', W_S));
fprintf(fid, sprintf('%0.2f Match Point W/P \n', HPW_match));
fprintf(fid, sprintf('%0.1f Reference Area (ft^2) \n', S));
fprintf(fid, sprintf('%0.1f Power (HP) \n\n', Power));

fprintf(fid, sprintf('%0.1f Loiter Velocity (Kts) \n', V_loiter));
fprintf(fid, sprintf('%0.3f Landing Mach Number \n', M_land));
fprintf(fid, sprintf('%0.3f Cruise Mach Number \n', M_cruise));
fprintf(fid, sprintf('%0.3f Loiter Mach Number \n', M_loiter));
fprintf(fid, sprintf('%0.3f Cruise CL \n', CL_cruise));
fprintf(fid, sprintf('%0.3f Loiter CL \n', CL_loiter));
fprintf(fid, sprintf('%0.0f Emergency Range (nmi) \n\n', Emergency_Dist));

fprintf(fid, sprintf('------------------------------------------------'));
fprintf(fid, sprintf('\nInput Parameters: \n \n'));

fprintf(fid, sprintf('%0.0f Takeoff Field Length (ft) \n',L_takeoff));
fprintf(fid,sprintf('%0.0f Landing Field Length (ft) \n \n',L_landing));

fprintf(fid, sprintf('%0.1f L/D at Cruise\n', LDC)); 
fprintf(fid, sprintf('%0.1f L/D at Loiter\n', LDL));
fprintf(fid, sprintf('%0.2f CL at Stall\n', Cl_stall)); 
fprintf(fid, sprintf('%0.2f CL Max TO\n', Clmax_to));
fprintf(fid, sprintf('%0.2f CL Max Land\n', Clmax_land));
fprintf(fid, sprintf('%0.2f CL Max Climb\n', Clmax_climb));
fprintf(fid, sprintf('%0.3f Zero-Lift Drag Coefficient\n', C_D0));
fprintf(fid, sprintf('%0.2f BSFC\n', bsfc));
fprintf(fid, sprintf('%0.2f Propeller Efficiency', eta));
fprintf(fid, sprintf('%0.1f Aspect Ratio\n', AR)); 
fprintf(fid, sprintf('%0.2f Oswald Efficiency Factor [e]\n', e)); 
fprintf(fid, sprintf('%0.3f Induced Drag Factor [k]\n', k));
fprintf(fid, sprintf('%0.3f Specific Heat Ratio [gamma] \n', gamma));
fprintf(fid, sprintf('%0.3f Gravity (ft/s^2) \n\n', g));

fprintf(fid, sprintf('%0.0f Range (Nm)\n', R));  
fprintf(fid, sprintf('%0.0f Stall Velocity (knots) \n', V_stall)); 
fprintf(fid, sprintf('%0.0f Cruise Velocity\n', V_cruise)); 
fprintf(fid, sprintf('%0.0f Approach Velocity\n', V_approach)); 
fprintf(fid, sprintf('%0.0f Climb Velocity\n', V_climb)); 
fprintf(fid, sprintf('%0.0f Rate of Climb (ft/min)\n', rate_climb)); 
fprintf(fid, sprintf('%0.0f Loiter Duration (min) \n', loiter_dur*60));
fprintf(fid, sprintf('%0.1f Max Load Factor (n max) \n ', n_max));
fprintf(fid, sprintf('%0.1f Min Load Factor (n min) \n ', n_min));
fprintf(fid, sprintf('%0.0f Number of Passengers\n ', passengers));
fprintf(fid, sprintf('%0.0f Number of Crew\n\n ', crew));

fprintf(fid, sprintf('%0.0f Cruise Altitude (ft) \n', altitude_ci)); 
fprintf(fid, sprintf('%0.0f Loiter Altitude (ft) \n', altitude_li)); 
fprintf(fid, sprintf('%0.0f Airfield Altitude (ft) \n', altitude_fi)); 
fprintf(fid, sprintf('%0.0f Climb Altitude (ft) \n\n', altitude_climbi));

fprintf(fid, sprintf('%0.1f V Tail Aspect Ratio(ft) \n', AR_v)); 
fprintf(fid, sprintf('%0.1f H Tail Aspect Ratio(ft) \n', AR_h)); 
fprintf(fid, sprintf('%0.1f Wing Taper Ratio \n', lambda)); 
fprintf(fid, sprintf('%0.1f V Tail Taper Ratio \n', lambda_v)); 
fprintf(fid, sprintf('%0.1f H Tail Taper Ratio \n', lambda_h)); 
fprintf(fid, sprintf('%0.1f Wing LE Sweep (deg) \n', sweep_LE));
fprintf(fid, sprintf('%0.1f V Tail LE Sweep (deg) \n', sweep_LE_v));
fprintf(fid, sprintf('%0.1f H Tail LE Sweep (deg) \n', sweep_LE_h));
fprintf(fid, sprintf('%0.2f V Tail Volume Coefficient (c_VT) \n', c_VT));
fprintf(fid, sprintf('%0.2f V Tail Volume Coefficient (C_HT) \n', c_HT));
fprintf(fid, sprintf('%0.0f Wing Dihedral (deg) \n', dihedral));
fprintf(fid, sprintf('%0.0f Fuselage Upsweep (deg) \n', upsweep));
fprintf(fid, sprintf('%0.2f Fuselage Diameter (ft) \n', D_fus));
fprintf(fid, sprintf('%0.2f Fuselage "a" Coefficient \n', a));
fprintf(fid, sprintf('%0.2f Fuselage "C" Coefficient \n', C));
fprintf(fid, sprintf('%0.2f Fuselage Fineness Ratio \n', fineness));

print(sprintf('%s/%s/Carpet Plot %s', folder_name, Trial_Name, ...
Trial_Name),'-dpng')

% fustabname = sprintf('%s/%s/Fuselage Sizing - %s.dat', folder_name, ...
%     Trial_Name,Trial_Name);
writetable(F,sprintf('%s/%s/Fuselage Sizing - %s.txt', folder_name, ...
    Trial_Name,Trial_Name),'WriteRowNames',true)
writetable(TSize,sprintf('%s/%s/Aircraft Sizing - %s.txt', folder_name, ...
    Trial_Name,Trial_Name),'WriteRowNames',true)
writetable(TBat,sprintf('%s/%s/Battery Specs - %s.txt', folder_name, ...
    Trial_Name,Trial_Name),'WriteRowNames',true)
end