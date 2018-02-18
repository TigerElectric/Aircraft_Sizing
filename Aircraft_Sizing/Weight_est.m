%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Leif Fredericks                                                       %%
%% AIAA 2017-2018 Hybrid-Electric General Aviation Aircraft (HEGAA)      %%
%% Whatever Bernardo is Calling This Directory                           %%
%% Feb. 15, 2018     
%% Dependencies: Ragone.m
%% Modified not yet
%% Aircraft Weight Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hardcoded variables to keep in mind
% N_prop for loiter is hardcoded to 0.7 as per Raymer, but may need to change
% 190 lb passengers and 30 lb for luggage

% MTOW[lb], R[nmi], V_cruise[Kts], V_stall[Kts], L_takeoff[ft],
% fpm_climb[ft/min], cruise_alt[ft], t_emergency_guess[Hrs]
function [W_TO, W_fuel, W_batt, W_empty, W_pay, W_batt_TO, W_batt_CM, W_batt_EM, t_EM] = ...
    Weight_est(MTOW, bsfc,AR, C_D0, Clmax_to, LDmax, N_prop, WP_cruise,...
    WP_climb, WP_TO, WP_match, WS, weight_max, R, V_cruise, V_stall, L_takeoff,...
    fpm_climb, cruise_alt, t_emergency_guess, passengers, crew, SHOW)
%% Old Guesses from Standalone (non-functional) version
% SHOW      =   0;      % SWITCH to show design points on non-log Ragone
% %%  STARRED ***** should be fed in as function from Prelim_Sizing
% % Performance Parameters
% MTOW        =   12000;   % FAR is 12,500 lb for FAR 23 classification*******
% bsfc        =   0.3;    % Aproximate based on Leif's engine # *************
% AR          =   7;      % Raymer guess ************************************
% C_D0        =   0.0275; % Estimate ****************************************
% Clmax_to    =   1.5;    % Good Est, NASA source [14]or[15] in Raymer ******
% LDmax       =   AR+10;  % early guess, should change **********************
% N_prop      =   0.8;    % Raymer, Reset for LOITER to 0.7 *****************
% 
% WP_cruise   =   17;     % rough guess, get from carpet plot ***************
% WP_climb    =   15;     % assume lower for now, get from carpet plot ******
% WP_TO       =   10;     % assume lower for now, get from carpet plot ******
% WS          =   25;     % based on 75 Kts stall, get from carpet plot *****
% weight_max  =   10^5;   % approximate maximum weight for iterations *******
% 
% % Design Requirements
% R           =   750;    % [nmi]     RFP ***********************************
% V_cruise    =   200;    % [Kts]     RFP ***********************************
% V_stall     =   61;     % [Kts]     FAR 23 ********************************
% L_takeoff   =   1800;   % [ft]      RFP ***********************************
% fpm_climb   =   1300;   % [ft/min]  RFP ***********************************
% cruise_alt  =   20000;  % [ft]      Design Book, see weight sources *******
% t_emergency_guess =   30 / 60;% [Hrs]     FAR-ish *******************************
% passengers  =   5;      %           RFP ***********************************
% crew        =   1;      %           RFP ***********************************
%%
% Derived Coefficients
V_climb     =   1.2*V_stall;                            % [Kts] Raymer 17.8.2
e           =   1.78*(1-0.045*AR^(0.68)) - 0.64;        % Raymer or Stengel
t_TO        =   L_takeoff/(V_climb*1.68781*0.7)/3600;   % rough takeoff time [hrs]
t_climb     =   cruise_alt / fpm_climb / 60;            % [Hrs]
WP          =   min([WP_cruise, WP_climb, WP_TO, WP_match]); % For Structural Considerations

% TAU = Vector of Discharge Times (vector can be any length)
TAU             =   [t_TO t_climb t_emergency_guess];         % [Hrs]
[PDens, EDens]  =   Ragone(TAU, SHOW);            % [W/kg] , [Whr/kg] 
PDens_HP        =   PDens*0.000608277;                  % [Hp/lb]
EDens_HP        =   EDens*0.000608277;                  % [Hp-hr/lb]
T = table(TAU', PDens, EDens); 
T.Properties.VariableNames = {'Time' 'WperKg' 'WhrperKg'};
% disp(T)
T2 = table(TAU', PDens_HP, EDens_HP);
T2.Properties.VariableNames = {'Time' 'HPperLB' 'HPHRperLB'};
% disp(T2)

%% Iterate Emergency Range Down
W_TO    =   weight_max; % Initial High Guess
while (W_TO > MTOW)

%% Fuel Weight Fraction
PI_beta     =   1;                % Running total BETA (Wi/W0)

% WARM-UP: Assumptions: 15 min
Beta_WU     = calculate_beta_hybrid('warmup', R, 15/60, V_cruise, V_stall, V_climb,...
    AR, e, C_D0, Clmax_to, bsfc, LDmax, N_prop, WP_cruise, PI_beta);
PI_beta     = Beta_WU;

% TAKE-OFF: 
Beta_TO     = calculate_beta_hybrid('takeoff', R, t_TO, V_cruise, V_stall, V_climb,...
    AR, e, C_D0, Clmax_to, bsfc, LDmax, N_prop, WP_cruise, PI_beta);
PI_beta     = Beta_TO*Beta_WU;

% CLIMB: Assumptions: Climb rate is constant
Beta_climb  = calculate_beta_hybrid('climb', R, t_climb, V_cruise, V_stall, V_climb,...
    AR, e, C_D0, Clmax_to, bsfc, LDmax, N_prop, WP_cruise, PI_beta);
PI_beta     = Beta_WU*Beta_TO*Beta_climb;

% CRUISE: 
Beta_cruise = calculate_beta_hybrid('cruise', R, 0, V_cruise, V_stall, V_climb,...
    AR, e, C_D0, Clmax_to, bsfc, LDmax, N_prop, WP_cruise, PI_beta);
PI_beta     = Beta_WU*Beta_TO*Beta_climb*Beta_cruise;

% LOITER: Assumptions: Propeller Efficiency is 0.7 as per Raymer
Beta_loiter = calculate_beta_hybrid('loiter', R, 0.75, V_cruise, V_stall, V_climb,...
    AR, e, C_D0, Clmax_to, bsfc, LDmax, 0.7, WP_cruise, PI_beta);
PI_beta     = Beta_WU*Beta_TO*Beta_climb*Beta_cruise*Beta_loiter;

% LAND: Assumptions: Approximately 5 min for landing
Beta_land   = calculate_beta_hybrid('land', R, 5/60, V_cruise, V_stall, V_climb,...
    AR, e, C_D0, Clmax_to, bsfc, LDmax, N_prop, WP_cruise, PI_beta);
PI_beta     = Beta_WU*Beta_TO*Beta_climb*Beta_cruise*Beta_loiter*Beta_land;

% FUEL WEIGHT FRACTION 
Wf_W0       =   1.06*(1-PI_beta);   % Weight Fraction of Fuel (Raymer 3.4.5)

%% Battery Weight Fraction

% EMERGENCY RANGE: Assumption: Cruise will be Most Efficient Speed by
% Engine Design, so Use Cruise Power
WBW0_EM         =   t_emergency_guess / (WP_cruise * EDens_HP(3));

% TAKEOFF POWER ASSIST:
if (WP_cruise > WP_TO)
    PTOB_W0     =   1/WP_TO - 1/WP_cruise;  %   Ratio of TO Power Supplied by Battery to Weight
    WBW0_TO     =   PTOB_W0 / PDens_HP(1);  
else
    disp('Take-Off Power is Less Than Cruise, so No Takeoff Battery')
    WBW0_TO     =   0;                      % Case where cruise is more limiting
end

% CLIMB POWER ASSIST:
if (WP_cruise > WP_climb)
    PCMB_W0     =   1/WP_climb - 1/WP_cruise;% Ratio of Climb Power Supplied by Battery to Weight
    WBW0_CM     =   PCMB_W0 / PDens_HP(2);  
else
    disp('Climb Power is Less Than Cruise, so No Climb Battery')
    WBW0_CM     =   0;                      % Case where cruise is more limiting
end

% CLIMB WEIGHT FRACTION
WB_W0   =   WBW0_EM + WBW0_TO + WBW0_CM;    % Total Weight Fraction of Batteries

%% Empty Weight Fraction
% From Raymer 6.3.3 for GA Single-Engine
a   =   -0.25;
b   =   1.18;
A   =   3.26;
c1  =   -0.20;
c2  =   0.08;
c3  =   0.05;
c4  =   -0.05;
c5  =   0.27; 

PHI     =   b * A^c2 * (1/WP)^c3 * WS^c4 * (V_cruise*1.68781)^c5; % V in [ft/s], lowest WP for structural limit

%% Iterative Solving

W_pay   =   (passengers+crew)*(190) + (passengers+crew)*30; % [lbs] Total Payload

w = linspace(1, weight_max, 1000);
raymer = @(W0) W0 - W0*(a + PHI*W0^c1) - W0*Wf_W0 - W0*WB_W0 - W_pay;

x0 = [1 weight_max];
try
    W_TO = fzero(raymer, x0);
catch
    errordlg('Aircraft weight not within specified range');
    error('Aircraft weight not within specified range');
end
t_emergency_final = t_emergency_guess;
t_emergency_guess = t_emergency_guess - 1/60;   % Take a minute off until becomes FAR 23
end
t_EM   =   t_emergency_final*60;
W_fuel              =   W_TO*Wf_W0;
W_batt              =   W_TO*WB_W0;
W_empty             =   W_TO*(a + PHI*W_TO^c1);
W_batt_EM           =   W_TO*WBW0_EM;
W_batt_TO           =   W_TO*WBW0_TO; 
W_batt_CM           =   W_TO*WBW0_CM;
if SHOW == 1
    figure()
    hold on;
    plot(w, w - w.*(a + PHI*w.^c1) - w.*Wf_W0 - w.*WB_W0 - W_pay);
    plot(W_TO,0,'r*');
    title('Raymer Equation')
    xlabel('Weight (lbs)');
    ylabel('N/A');
    hold off;
end
end






