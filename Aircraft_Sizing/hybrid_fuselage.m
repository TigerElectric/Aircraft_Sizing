% Upsweep in [deg], Weights in [lbs], measurements in [ft]
function [L_fus, L_HT, L_VT, V_fus, D_fus_rec, S_wet_fus, L_fus_angled, V_batt] = ...
    hybrid_fuselage(a, W_TO, W_batt, C, fineness, L_nose, upsweep, D_fus)
% Going to Assume average 250 lb/ft3 for battery
V_batt = W_batt / 250; % Volume of batteries alone
V_batt_all = V_batt*2; % Include twice the volume for insulation, siding, etc
L_fus = a*(W_TO-W_batt)^C; % ft
D_fus_rec = L_fus/fineness; % recommended fuselage diamter ft
L_fus = L_fus + V_batt_all/(pi*(D_fus/2)^2);
L_HT = 0.60 * L_fus; % ft, H tail moment arm : aft engine per Raymer 160]
L_VT = L_HT; 

% All the following is Raymer 7.9-7.10
% Top View Calculation
A_top = D_fus*(L_fus-L_nose)+(.5*(D_fus + 0.2*D_fus)*L_nose); 
                     % ft, trapezoid assuming that the diameter of the prop cone is 1/5 
                     % diameter of fuselage from top view
% Side View Calculation
L_fus_angled = D_fus*0.75/ tand(upsweep); % ft, horizontal distance from rear
                                       % where fuselage angles upward
                                       % ASSUMING tail terminates 1/4 of
                                       % the diameter above the centerline
A_side = D_fus*(L_fus-L_fus_angled-L_nose)+(.5*D_fus*L_fus_angled)...
    +(.5*(D_fus + 0.5*D_fus)*L_nose);
                     % ft, assuming that the diameter of the prop cone is 1/2 
                     % diameter of fuselage from side view                   
% Volume Calculation
V_fus = 3.4 * A_top * A_side / (4*L_fus);

% Wetted Area Calculation
S_wet_fus = 3.4*(A_top * A_side)/2;
end