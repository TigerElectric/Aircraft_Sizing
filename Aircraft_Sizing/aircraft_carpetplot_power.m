%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bernardo Pacini                                                       %%
%% Hybrid-Electric General Aviation Aircraft                             %%
%% Preliminary Design Calculations - Aircraft Carpet Plot                %%
%% Created: 1/23/2018                                                    %%
%% Modified: 2/15/2018                                                   %%
%% Dependencies:                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [HPW_takeoff_ws, HPW_climb_ws, HPW_cruise_ws, HPW_turning_max_ws, ...
    HPW_turning_min_ws, HPW_endurance_ws, wing_loading] = aircraft_carpetplot_power(V_cruise, AR, eta, WP_guess, WS_guess,...
    altitude_ci, altitude_fi, altitude_climbi, V_approach, ...
    V_stall, Clmax_to, Clmax_climb, Clmax_land, L_takeoff, L_landing, V_climb, ...
    rate_climb, theta_app, C_D0, gamma, g, ...
    carpet_x_lim, carpet_y_lim, LD, LDC, n_max, n_min, k, e)

altitude_c = altitude_ci*0.3048;
altitude_climb = altitude_climbi*0.3048;
altitude_f = altitude_fi*0.3048;
[airDens_c, airPres_c, temp_c, soundSpeed_c] = Atmos(altitude_c); % SI
[airDens_f, airPres_f, temp_f, soundSpeed_f] = Atmos(altitude_f);
[airDens_climb, airPres_climb, temp_climb, soundSpeed_climb] = ...
    Atmos(altitude_climb);
[airDens_sl, airPres_sl, temp_sl, soundSpeed_sl] = Atmos(0);

% Convert values from SI to Imperial
[airDens_ci, airPres_ci, temp_ci, soundSpeed_ci] = ...
    convert_to_imperial(airDens_c, airPres_c, temp_c, soundSpeed_c);
[airDens_fi, airPres_fi, temp_fi, soundSpeed_fi] = ...
    convert_to_imperial(airDens_f, airPres_f, temp_f, soundSpeed_f);
[airDens_climbi, airPres_climbi, temp_climbi, soundSpeed_climbi] = ...
    convert_to_imperial(airDens_climb, airPres_climb, temp_climb, ...
    soundSpeed_climb);
[airDens_sli, airPres_sli, temp_sli, soundSpeed_sli] = ...
    convert_to_imperial(airDens_sl, airPres_sl, temp_sl, soundSpeed_sl);

sigma = airDens_fi/airDens_sli;
V_stall = V_stall * 1.68781; %convert to ft/s
dHdt = rate_climb/60; % ft/s
V_climb = V_climb*1.68781; % ft/s

e = 1.78*(1-0.045*AR^(0.68)) - 0.64;
k = 1/(pi*AR*e); % induced drag correction factor, 


figure()
hax=axes; 
title('Constraint Plane (T/W - W/S)');
xlabel('Wing Loading [W_g/S], lb/ft^2');
ylabel('Thrust Loading [T_0/W_g]');
hold on;
% Plotted later for cosmetic reasons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stall DONE
WS_stall = ((V_stall^2)*airDens_sli*Clmax_land)/(2*g);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take-off Parameter DONE
% (W/S) = (TOP)*sigma*Cl_TO*(hp/W)
TOP = 250;
WS = linspace(1,200);
HPW_takeoff = @(ws) (TOP * sigma * Clmax_to)./ws;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIELD LENGTHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mum = 0.5
% S_A = 1800;
% kt = @(mu) (550*eta)/(1.2*V_cruise*.7*WP_guess) - mu
% ka = @(mu) ((airDens_fi/g)/(2*WS_guess))*(mu*Clmax_to - C_D0 - k*Clmax_to^2)
% disp(kt(mum))
% disp(ka(mum))
% s_g =  (1/(2*g*ka(mum))) * log((kt(mum) + ka(mum)*(1.2*V_stall))/(kt(mum) + ka(mum)*0))
% %test = s_g(0.03)


%%%%%%%%%%%%%%%%%%%%%%%%q%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Climb-performance DONE
gamma = asind(dHdt/V_climb);
HPW_climb = @(ws) (550*eta)*(1/(sind(gamma) + (C_D0/(.5*Clmax_to))))...
    *sqrt((airDens_ci*Clmax_climb)./(2*g.*ws));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cruise-performance DONE
V_cruise = V_cruise * 1.68781; %convert to ft/s
HPW_cruise = @(ws) 550./(0.5*(airDens_ci/g)*V_cruise^3*C_D0./ws ...
     + k /(0.5*(airDens_ci/g)*V_cruise) .* ws);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Turning Flight DONE
q_c = dynamic_pressure(airDens_ci, V_cruise, g);
HPW_turning_max = @(ws) (550*eta)./(V_cruise*(((n_max^2 * k * ws)./(q_c))...
    + ((q_c*C_D0)./(ws))));
HPW_turning_min = @(ws) (550*eta)./(V_cruise*(((n_min^2 * k * ws)./(q_c))...
    + ((q_c*C_D0)./(ws))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Endurance DONE
HPW_endurance = @(ws) ((3*C_D0/k)^(3/4))*(sqrt(airDens_ci)*550)./...
    (sqrt(2*ws)*4*C_D0*sqrt(g));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Landing DONE
WS_landing = ((sqrt(L_landing/0.265)*1.687)^2)*...
    (airDens_sli/g*Clmax_land/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot everything on carpet plot
% plot rectangles first so they're in back
rectangle('Position', [WS_stall carpet_y_lim(1) ...
    carpet_x_lim(2) carpet_y_lim(2)], 'FaceColor', [1 1 0]);
rectangle('Position', [WS_landing carpet_y_lim(1) ...
   (carpet_x_lim(2)) carpet_y_lim(2)], 'FaceColor', [0 1 1]);
% Stall
line([WS_stall WS_stall],get(hax,'YLim'),'Color',[1 1 0]);
% Landing
line([WS_landing WS_landing], get(hax,'YLim'),'Color',[0 1 1]);
alpha(0.5)

% Takeoff
area(WS, HPW_takeoff(WS), carpet_y_lim(2), 'FaceColor', 'b');
alpha(0.5); % transparency
% Climb
area(WS, HPW_climb(WS), carpet_y_lim(2), 'FaceColor', 'r');
alpha(0.5); % transparency
% Cruise
area(WS, HPW_cruise(WS), carpet_y_lim(2), 'FaceColor', 'g');
alpha(0.5); % transparency
% Endurance
area(WS, HPW_endurance(WS), carpet_y_lim(2), 'FaceColor', 'm');
alpha(0.5); % transparency
% Load Factors
plot(WS, HPW_turning_max(WS), 'r-.');
alpha(0.5); % transparency
plot(WS, HPW_turning_min(WS), 'r--');
alpha(0.5); % transparency
title('Constraint Plane (T/W - W/S)');
xlabel('Wing Loading [W_g/S], lb/ft^2');
ylabel('Thrust Loading [T_0/W_g]');
xlim(carpet_x_lim)
ylim(carpet_y_lim)

% Takeof
plot(WS, HPW_takeoff(WS), 'b');
% Climb
plot(WS, HPW_climb(WS), 'r');
% Cruise
plot(WS, HPW_cruise(WS), 'g');
% Loiter
plot(WS, HPW_endurance(WS), 'm');
alpha(0.5); % transparency

title('Constraint Plane (W/P - W/S)');
xlabel('Wing Loading [W_g/S], lb/ft^2');
ylabel('Loading Power [W_g/P_{hp}]');
legend('Stall', 'Landing','Takeoff','Climb','Cruise','Endurance',...
    'n_{max}','n_{min}');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% THRUST AND SURFACE AREA VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prompt = {'Power Loading [W_g/P]:','Wing Loading [W_g/S]:'};
dlg_title = 'User Input';
num_lines = 1;
defaultans = {'Value...','Value...'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

power_loading = str2num(answer{1});
wing_loading = str2num(answer{2});
plot(wing_loading,power_loading,'b.', 'MarkerSize', 20);

HPW_takeoff_ws     = HPW_takeoff(wing_loading);
HPW_climb_ws       = HPW_climb(wing_loading);
HPW_cruise_ws      = HPW_cruise(wing_loading);
HPW_turning_max_ws = HPW_turning_max(wing_loading);   
HPW_turning_min_ws = HPW_turning_min(wing_loading);
HPW_endurance_ws   = HPW_endurance(wing_loading);
end