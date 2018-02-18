%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Leif Fredericks                                                       %%
%% AIAA 2017-2018 Hybrid-Electric General Aviation Aircraft (HEGAA)      %%
%% Whatever Bernardo is Calling This Directory                           %%
%% Dec. 3 2017      
%% Dependencies: Ragone.m | 
%% Modified not yet
%% Tests Modular Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Test Ragone.m 
SHOWRAGONE      =   0;      % SWITCH to show design points on non-log Ragone
t1              =   60;     % [s]
t2              =   5;      % [min]
t3              =   1;      % [hr]
% Vector of as many timescales as you want to investigate
TAU             =   [t1/3600 t2/60 t3];      % Battery Discharge Times (Hrs)
[PDens, EDens]  =   Ragone(TAU, SHOWRAGONE); % [W/Kg ; Whr/Kg]
T = table([t1; t2; t3], PDens, EDens); 
T.Properties.VariableNames = {'Time' 'WperKg' 'WhrperKg'};
disp(T)
T2 = table([t1; t2; t3], PDens*0.000608277, EDens*0.000608277);
T2.Properties.VariableNames = {'Time' 'HPperLB' 'HPHRperLB'};
disp(T2)
%% Test Beta with 6 person Variant
% % R (nmi), E (hrs), V (knots), bsfc (lbfuel/hr/Hp), Climb_time (hrs)
% function beta = calculate_beta_hybrid(ID, R, E, V_cruise, V_stall, V_climb,...
%     AR, e, C_D0, Clmax_to, bsfc, LDmax, N_prop, WP_cruise)
bsfc = 0.3;
R = 750; % nmi
V_cruise = 200; % knots
V_stall     =   61; % knots
V_climb     = 1.2*V_stall; % knots
AR = 7; % guess
e = 1.78*(1-0.045*AR^(0.68)) - 0.64;
C_D0 = 0.0275; % est
Clmax_to = 1.5; % good est
LDmax = AR+10; % early guess
N_prop = 0.8; % Reset for Loiter to 0.7
WP_cruise = 12; % rough guess
TO_time = 1800 / (V_climb*1.68781*0.7) /3600; % rough takeoff time [hrs]
fpm_climb = 1300 ; % feet/minute climb requirement
cruise_alt = 20000; % feet
climb_time = cruise_alt / fpm_climb / 60; % hours


% Warmup
PI_beta = 1;

Beta_WU = calculate_beta_hybrid('warmup', R, 15/60, V_cruise, V_stall, V_climb,...
    AR, e, C_D0, Clmax_to, bsfc, LDmax, N_prop, WP_cruise, PI_beta);
PI_beta = Beta_WU;

Beta_TO = calculate_beta_hybrid('takeoff', R, TO_time, V_cruise, V_stall, V_climb,...
    AR, e, C_D0, Clmax_to, bsfc, LDmax, N_prop, WP_cruise, PI_beta);
PI_beta = Beta_TO*Beta_WU;

Beta_climb = calculate_beta_hybrid('climb', R, climb_time, V_cruise, V_stall, V_climb,...
    AR, e, C_D0, Clmax_to, bsfc, LDmax, N_prop, WP_cruise, PI_beta);
PI_beta = Beta_WU*Beta_TO*Beta_climb;

Beta_cruise = calculate_beta_hybrid('cruise', R, 0, V_cruise, V_stall, V_climb,...
    AR, e, C_D0, Clmax_to, bsfc, LDmax, N_prop, WP_cruise, PI_beta);
PI_beta = Beta_WU*Beta_TO*Beta_climb*Beta_cruise;

Beta_loiter = calculate_beta_hybrid('loiter', R, 0.75, V_cruise, V_stall, V_climb,...
    AR, e, C_D0, Clmax_to, bsfc, LDmax, 0.7, WP_cruise, PI_beta);
PI_beta = Beta_WU*Beta_TO*Beta_climb*Beta_cruise*Beta_loiter;

Beta_land = calculate_beta_hybrid('land', R, 5/60, V_cruise, V_stall, V_climb,...
    AR, e, C_D0, Clmax_to, bsfc, LDmax, N_prop, WP_cruise, PI_beta);
PI_beta = Beta_WU*Beta_TO*Beta_climb*Beta_cruise*Beta_loiter*Beta_land;

%% Hardcoded variables to keep in mind

SHOW      =   0;      % SWITCH to show design points on non-log Ragone
%%  STARRED ***** should be fed in as function from Prelim_Sizing
% Performance Parameters
MTOW        =   7000;   % FAR is 12,500 lb for FAR 23 classification*******
bsfc        =   0.3;    % Aproximate based on Leif's engine # *************
AR          =   7;      % Raymer guess ************************************
C_D0        =   0.0275; % Estimate ****************************************
Clmax_to    =   1.5;    % Good Est, NASA source [14]or[15] in Raymer ******
LDmax       =   AR+10;  % early guess, should change **********************
N_prop      =   0.8;    % Raymer, Reset for LOITER to 0.7 *****************

WP_cruise   =   17;     % rough guess, get from carpet plot ***************
WP_climb    =   13;     % assume lower for now, get from carpet plot ******
WP_TO       =   14;     % assume lower for now, get from carpet plot ******
WP_match    =   15;     % will
WS          =   25;     % based on 75 Kts stall, get from carpet plot *****
weight_max  =   10^6;   % approximate maximum weight for iterations *******

% Design Requirements
R           =   750;    % [nmi]     RFP ***********************************
V_cruise    =   200;    % [Kts]     RFP ***********************************
V_stall     =   61;     % [Kts]     FAR 23 ********************************
L_takeoff   =   1800;   % [ft]      RFP ***********************************
fpm_climb   =   1300;   % [ft/min]  RFP ***********************************
cruise_alt  =   20000;  % [ft]      Design Book, see weight sources *******
t_emergency_guess =   60 / 60;% [Hrs]     FAR-ish *******************************
passengers  =   5;      %           RFP ***********************************
crew        =   1;      %           RFP ***********************************


% N_prop for loiter is hardcoded to 0.7 as per Raymer, but may need to change
% 190 lb passengers and 30 lb for luggage

% MTOW[lb], R[nmi], V_cruise[Kts], V_stall[Kts], L_takeoff[ft],
% fpm_climb[ft/min], cruise_alt[ft], t_emergency_guess[Hrs]


[W_TO, W_fuel, W_batt, W_empty, W_pay, W_batt_TO, W_batt_CM, W_batt_EM, t_EM]...
    = Weight_est(MTOW, bsfc,...
AR, C_D0, Clmax_to, LDmax, N_prop, WP_cruise, WP_climb, WP_TO, WP_match, WS,...
weight_max, R, V_cruise, V_stall, L_takeoff, fpm_climb, cruise_alt,...
t_emergency_guess, passengers, crew, SHOW)

t_EM/60*V_cruise




