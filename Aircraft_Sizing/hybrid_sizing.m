
function [TSize] = hybrid_sizing(S, AR, lambda, sweep_LE, AR_v, AR_h, lambda_v, ...
    lambda_h, sweep_LE_v, sweep_LE_h, c_VT, c_HT, L_VT, dihedral,...
    L_fus, D_fus) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Size the Main Wing
[b, C_root, C_tip, C_bar, Y_bar] = ...
    wing_sizing(S, AR, lambda); 
% Find Main Wing Sweeps
[sweep_25, sweep_TE] = ...
    sweep(b, sweep_LE, C_root, C_tip);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Size Vertical Tail, Assuming Coventional Arrangement
S_VT        =   c_VT * b * S / L_VT;    % Raymer 159 Vertical Tail Area
[h_v, C_root_v, C_tip_v, C_bar_v, Y_bar_v] = ...
    wing_sizing(S_VT, AR_v, lambda_v);  % Set Geometric Parameters for VT
dihedral_v  =   90;                     % Dihedral of a Vertical Tail    
% Vertical Tail Sweeps, recal h_v is 1/2 the "span"
[sweep_25_v, sweep_TE_v] = ...
    sweep(2*h_v, sweep_LE_v, C_root_v, C_tip_v);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Size Horizontal Tail, Assuming Coventional Arrangement
L_HT        =   L_VT;   % Conservative cause HT should be slightly back from VT
S_HT        =   c_HT * C_bar * S / L_HT;    % Raymer 159 Horizontal Tail Area
[b_h, C_root_h, C_tip_h, C_bar_h, Y_bar_h] = ...
     wing_sizing(S_HT, AR_h, lambda_h);  % Set Geometric Parameters for HT
dihedral_h  =   dihedral;               % Assume same as wing    
% Horizontal Tail Sweeps
[sweep_25_h, sweep_TE_h] = ...
    sweep(b_h, sweep_LE_h, C_root_h, C_tip_h);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Position and Length Parameters
% Vertical Tail Position
X_v         = L_fus - C_root_v - 1;     % Bring it 1 foot foreward from back for now
Z_v         = D_fus/4;                  % Same est as fuselage that aft 
                                        % terminates 1/2 radius avove
X_v_AC      = X_v + Y_bar_v*tand(sweep_LE_v) + C_bar_v/4;
Z_v_AC      = Z_v + Y_bar_v;
% Wing Position
X_w_AC      = X_v_AC - L_VT;            % 1/4 MAC of wing, ie x aero center
X_w         = X_w_AC - Y_bar*tand(sweep_LE) - C_bar/4;
Z_w         = -D_fus/2 +.5;
Z_w_AC      = Z_w + Y_bar * tand(dihedral);
% Horizontal Tail Position
X_h         = L_fus - C_root_h;         % Back right to aft
Z_h         = D_fus/4;                  % Same est as fuselage that aft 
                                        % terminates 1/2 radius avove
                                        % centerline
X_h_AC      = X_h + Y_bar_h*tand(sweep_LE_h) + C_bar_h/4;
Z_h_AC      = Z_h + Y_bar_h*tand(dihedral_h);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Parameters = {'Span'; 'Surface Area'; 'AR'; 'Root Chord'; 'Tip Chord'; ...
    'Taper Ratio'; 'LE Sweep'; 'Dihedral'; 'X Position'; 'Z Position';...
    'X Aerodynamic Center';'Z Aerodynamic Center';'Mean Aerodynamic Chord';...
    'Spanwise MAC Location';'TE Sweep';'Quarter Chord Sweep'};
Wing = [b; S; AR; C_root; C_tip; lambda; sweep_LE; dihedral; X_w; Z_w;...
    X_w_AC; Z_w_AC;C_bar;Y_bar;sweep_TE;sweep_25];
Horizontal_Tail = [b_h; S_HT; AR_h; C_root_h; C_tip_h; lambda_h;...
    sweep_LE_h; dihedral_h; X_h; Z_h; X_h_AC; Z_h_AC;C_bar_h; Y_bar_h;...
    sweep_TE_h;sweep_25_h];
Vertical_Tail = [h_v; S_VT; AR_v; C_root_v;C_tip_v; lambda_v; ...
    sweep_LE_v; dihedral_v; X_v; Z_v; X_v_AC; Z_v_AC;C_bar_v;Y_bar_v;...
    sweep_TE_v;sweep_25_v];
T = table(Wing, Horizontal_Tail, Vertical_Tail, 'RowNames', Parameters);
% disp(T);
TSize = T;
end
% Size wing parameters based on Raymer 192 assuming trapezoidal wing
% Works for VT height cause of def of AR, so no need to half anything
function [b, C_root, C_tip, C_bar, Y_bar] = ...
    wing_sizing(S, AR, lambda)
b = sqrt(S * AR);
C_root = 2 * S / (b *(1+lambda));
C_tip = lambda * C_root;
C_bar = (2/3) * C_root * (1 + lambda + lambda^2) ...
    / (1 + lambda);
Y_bar = (b / 6) * (1 + 2*lambda) / (1 + lambda);
end
% Determine Trailing edge and Quarter Chord Sweeps Geometrically
function [sweep_c, sweep_TE] = sweep(b, sweep_LE, C_root, C_tip)
z = b/2 * tand(sweep_LE);
x = z - (C_root / 4);
y = x + (C_tip / 4);
sweep_c = atand(y/(b/2));
sweep_TE = atand((z + C_tip - C_root) / (b/2));
end
