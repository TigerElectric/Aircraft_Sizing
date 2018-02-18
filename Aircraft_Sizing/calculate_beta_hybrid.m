%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bernardo Pacini & Nathan Wei                                          %%
%% MAE 332 - Aircraft Design                                             %%
%% Preliminary Design Calculations                                       %%
%% 7 March 2017  
%% Dependencies: none
%% Modified 07/02/2017 by Leif Fredericks
%% Calculate Beta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% R (nmi), E (hrs), V (knots), bsfc (lbfuel/hr/Hp), Climb_time (hrs)
function beta = calculate_beta_hybrid(ID, R, E, V_cruise, V_stall, V_climb,...
    AR, e, C_D0, Clmax_to, bsfc, LDmax, N_prop, WP_cruise, PI_beta)

LDC                 =   LDmax; % Raymer 41, cruise is at max L/D
LDL                 =   LDmax * 0.866; % Raymer 41, loiter is at 86.6% of max L/D
V_prop_conversion =   V_cruise * 1.68781; % Feet / second from knots
C_cruise  =   bsfc * V_prop_conversion / 550 / N_prop; % Raymer 3.4.3, [lb/(lb-hr)] <- thrust specific
TW_cruise          =  550*N_prop/V_prop_conversion / WP_cruise;

switch ID
    
    case 'warmup'
        WP_idle =   WP_cruise / 0.2; % Diesel (conservative) estimate assumes 600rpm idle
        WU      =   E; % Time spent warming up [hrs]
        beta    =   1 - (1/WP_idle * bsfc * WU); % [hp/lb][lb/(hp*hr)][hr]
    
    case 'takeoff'
        TO          = E; % Time spent lifting off [hrs]
        % Lots from http://www.fzt.haw-hamburg.de/pers/Scholz/HOOU/AircraftDesign_5_PreliminarySizing.pdf
        e           = 0.7; % Reduces for flaps deployed 
        V_2         = 1.2*V_stall; % Final Velocity of take-off segment (hamburg)
        CL          = Clmax_to; % * (V_stall/V_2)^2;
        dCD_flap    = 0.05*CL - 0.055; %  http://www.fzt.haw-hamburg.de source
        dCD_gear    = 0.015; % http://www.fzt.haw-hamburg.de
        
        V_prop_conversion =   V_2 * 1.68781; % Feet / second from knots
        C  =   bsfc * V_prop_conversion / 550 / N_prop; % Raymer 3.4.3, [lb/(lb-hr)] <- thrust specific
        LDTO        = CL / ((C_D0+dCD_flap+dCD_gear + ...
                        (CL^2 / (3.14*AR*e))));
        
        beta = 1/(exp((TO*C)/LDTO));
    case 'climb' % assume climb clean, no flaps
        CLIMB       = E; % Time spent climbing [hrs]
        V_climb     = V_climb * 1.68781; % [FPS] from knots
        K           = 1/(pi*AR*e); 
        TW          =  550*N_prop/V_climb / WP_cruise; % to T/W, can use WP cruise cause what engine sees
        CL          = 1/(2*K) *(-TW + sqrt(TW^2+12*C_D0*K));
        C           =   bsfc * V_climb / 550 / N_prop; % Raymer 3.4.3, [lb/(lb-hr)] <- thrust specific      
        LDCL        = CL / ((C_D0 + (CL^2 / (3.14*AR*e))));
        beta = 1/(exp((CLIMB*C)/LDCL));
        
    case 'cruise'
        % Breguet Range Equation
        % R = (V/C) * (L_D) * ln(Wi/Wf) %lbfuel/h/lbt
 
        beta = 1/(exp(R*((C_cruise/V_cruise))/(LDC)));
        
    case 'loiter'
        % Use cruise speed for conservative estimate 
        
        C  =   bsfc * V_prop_conversion / 550 / 0.7; % drop N_prop to 0.7 (Raymer)
        beta = 1/(exp((E*C)/LDL));
    
    case 'descent'
        WP_idle =   WP_cruise / 0.2; % Diesel (conservative) estimate assumes 600rpm idle
        DESC      =   E; % Time spent landing [hrs]
        beta    =   1 - (bsfc * DESC)/(WP_idle * PI_beta); % [hp/lb][lb/(hp*hr)][hr]
        
    case 'land'
%         beta = 0.9725; % estimate
        WP_idle =   WP_cruise / 0.2; % Diesel (conservative) estimate assumes 600rpm idle
        WP_idle = WP_idle/2; % Safety Factor of 2 on power in case need more for landing
        APP      =   E; % Time spent landing [hrs]
        beta    =   1 - (bsfc * APP)/(WP_idle * PI_beta); % [hp/lb][lb/(hp*hr)][hr]
        
end

end