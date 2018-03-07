%% Depends on Atmos.m
clear all;
% close all;

%% Output Controls
plotPV  =   1; % Plots PV Diagram of Cycle
Table   =   1; % Plots Table of P,V,T at All States
%% Design Decision
Diesel  =   1; % Set to 1 if using diesel CYCLE (with Jet A-1 or Diesel), 0 for avgas
%% Most Relevant Design Factors (NEED TO ENSURE NetTurbo > 0!!!)
%% Input Block
Power   =   732;        % Horsepower
Alt     =   6096;       % Cruise Altitude (20,000 ft = 6096) [m]
Boost   =   4.8;        % Ratio of Intake to Atm Pressure(Pi/Patm) 1.1254
if(Alt == 0)Boost = 0.43*Boost;end % Use WASTEGATE to only have boost at alt
if(Boost < 1) Boost = 1; end; % just in case ratio drops it too low
xq      =   0.8;        % Fraction of Heat Input that is ISOCHORIC
Eps     =   20;          % Effective Compression Ratio [] i like 20
Sigma   =   1.2;          % 1-1' Additional Expansion Ratio [] 1.4
Vd      =   16.9;          % Displacement volume (Liters)
Lean    =   0.89;       % Fraction of stoich fuel to air ratio, .5 gud, 0.606 max for diesel

Vd      =   Vd * 0.001; % Covcert to [m^3]
Plim    =   1.034E7;    % Maximum Pressure (1.034E7 for avgas ; 2E7 diesel) [Pa] 
LHV     =   4.354E7;    % Lower Heating Value of Avgas [J/Kg]
Stoich  =   0.0673;      % Avgas approximate stoichiometric ratio 1:15
redline =   5500;
if (Diesel == 1) 
    LHV = 42.8E6; 
    xq = 0.2; 
    redline = 3000;
    Plim = 2E7;
    Lean    =   0.606;
end % Assume appx 25% isochoric in modern CI

Tlim    =   4000;       % Maximum Temperature [K] (Probably 2600K)


%% Other Design Constants

PSL             =   101300;     % Seal Level Atmospheric Pressure
R               =   287;        % Ideal Gas Constant [J/kg-K]
y               =   1.4;        % Adiabatic Index
pi              =   3.1415;     % Pi
N_c             =   0.70;       % 70% Compresor Efficiency, guess
N_t             =   0.80;       % 80% Turbine Efficiency, guess 
N_m             =   0.95;       % 95% Mechanical Efficiency, guess, not used here, but negligible(mention)
N_t             =   N_t*N_m;    % Account for Mechanical Efficiency
theory2actual   =   0.80;       % Liberal side of Mueller range(.7), several sources(.8)
N_power         =   0.90;       % Brake vs indicated (UofIdaho Slides)
N_comb          =   0.95;       % combustion efficiency http://slideplayer.com/slide/10276479/# slide 49

%% Derived Design Factors, DO NOT MODIFY BELOW THIS POINT

[airDens,Patm,Tamb,soundSpeed] = Atmos(Alt); % [kg/m^3, N/m^2, K, m/s]
Pi          =   Patm * Boost; % Intake Pressure [Pa]
Pback       =   2/3 * (Pi - Patm); % Based on 2:1 back:boost at SL (cit need) [Pa]
Pe          =   Patm + Pback; % Absolute Exhaust Pressure [Pa]
B           =   Pi / Pe; % Intake to Exhaust Pressure Ratio []
rc          =   Eps * Sigma; % Full Expansion/Compression Ratio []
rp          =   Pi/Patm; % Pressure ratio for turbocharger (shuld equal Boost)
n           =   y; % Polytropic Exponent Assuming Isentropic 
Cv          =   R / (y-1); % Isochoric Heat Capacity [J/kg-K]
Cp          =   R * y / (y-1); % Isobaric Heat Capacity [J/kg-K]
FA          =   Stoich * Lean;  % Fuel/Air Ratio, Best Will Likely Be in Range (0.055-0.065) 
LHV_kwhg    =   LHV * 2.77778e-10; % Lower Heating Value in [kWh/g]

%% State 8 Start Intake Stroke
% Assume can find intake temperature from isentropic compression from
% atmosphere through compressor to intake pressure
V8  =   Vd / (rc - 1); % Volume at TDC [m^3]
[T8, P8, V8, wc] = JANAF(Tamb, Patm, 1, V8, 0, 0, 0, 5, rp, N_c, N_t, 0);

%% State 1 Effective End Intake Stroke / Start Compression Stroke

% Assume this is an open system and the intake temperature is the same
% during intake stroke
P1  =   Pi; % Intake is Isobaric [Pa]
V1  =   V8 * Eps; % Definition of Epsilon [m^3]
T1  =   T8; % Open System So Intake Manifold Temp [K]
W81 =   Pi * (V1 - V8);
m   =   P1*V1/(R*T1); % Mass of Air After Intake

%% State 2 End Compression Stroke / Start Isochoric Heat 

[T2, P2, V2, W12] = JANAF(T1, P1, m, V1, 0, Eps, 0, 1);

%% State 3 End Isochoric Heat / Start Power Stroke and Isobaric Heat

Qin =   FA*LHV*m*N_comb;
Qin23 = Qin*xq;
[T3, P3, V3, W23] = JANAF(T2, P2, m, V2, Qin23, 0, 0, 2);

%% State 4 End Isobaric Heat / Continue Power Stroke

Qin34 = Qin*(1-xq);
[T4, P4, V4, W34] = JANAF(T3, P3, m, V3, Qin34, 0, 0, 3);

%% State 5 End Power Stroke / Start Isochoric Heat Rejection

BetaBar = V4/V3;
re = Eps*Sigma/BetaBar;
[T5, P5, V5, W45] = JANAF(T4, P4, m, V4, 0, 0, re, 4);
if (P5 < Pe)
    error('Engine Tries to Expand below Exhaust Pressure')
end

%% State 6 End Isochoric Heat Rejection / Start Exhaust Stroke

V6  =   V5; % Isochoric Process [m^3]
P6  =   Pe; % Exhaust at Back Pressure [Pa]
T6  =   T5 * (Pe / P5); % Ideal Gas Law

%% State 6' (61) Available Expansion

V61         =   V5 * ((P5/Pe)^(1/y)); % Isentropic Relation [m^3]
P61         =   Pe; % Exhaust at Back Pressure [Pa]
T61         =   T5 * ((V5/V61)^(y-1)); % Isentropic Relation [K]

%% State 6'' (62) Actual Turbine Exit
% Use wc, and recall it is supposed to be mass-normalized
[T62, P62, V62, wt] = JANAF(T5, P5, m, V5, 0, 0, 0, 6, 0, N_c, N_t, wc);

% T62         =   T5 - Tamb/N_c*((rp)^((n-1)/n)-1);
% P62         =   P5*(1 - 1/N_c/N_t * Tamb/T5 *((rp)^((n-1)/n)-1))^(n/(n-1));

% Purely a check. These should all be equal for turbocharger to be valid.
% ISENTROPIC
WT      =   N_t*m*Cp*T5*(1-(P62/P5)^((n-1)/n));
WT2     =   -m*Cp*(T62-T5);
WC      =   m*Cp*(T8-Tamb);
WC2     =   m*Cp*Tamb/N_c * (rp^((n-1)/n)-1);


%% State 7 End Exhaust Stroke / Start Isochoric Gas Exchange

P7  =   Pe; % Exhaust at Back Pressure [Pa]
V7  =   V2; % TDC [m^3]
T7  =   T6; % Open System, doesn't really matter [K]
W67 =   Pe * (V7 - V6);

%% Cycle Properties

W_ideal     =   W12 + W34 + W45 + W67 + W81; % Net Work of System
W           =   W_ideal*N_power; % Brake Work
Eta_ideal   =   W/Qin; % Ideal Cycle Efficiency, Accounts for Brake vs Indicated already 
Eta         =   theory2actual * Eta_ideal; 
Eff         =   Eta*100;
imep        =   W/Vd; % Indicated mean effective pressure (Pa)
bmep        =   N_power * imep; % Brake mean effective pressure
bsfc_met    =   1/(Eta*LHV_kwhg); % Metric Brake Specific Fuel Consumption [g/kWh]
bsfc_bhp    =   0.001644 * bsfc_met; % BSFC [lb / (hp*h)]
torque      =   bmep * Vd / (4*pi); % Brake Torque [N-m]
P_watts     =   Power * 746; % Power Specified in Watts
RPM         =   60* P_watts / (2*pi*torque); % want 5-7 K based http://slideplayer.com/slide/10276479/#
MaxPower    =   (2*pi*torque) * redline / 60 / 746;
%%

if (plotPV ==1)
figure
hold on
plot(V8,P8,'r*')
plot([V8 V1], [P8 P1])
plot(V1,P1,'r*')
x12 = linspace(V1,V2,100);
P2graph = [];
for i = 1:100
[T2graph, P2graphthis, V2graph] = JANAF(T1, P1, m, V1, 0, V1/x12(i), 0, 1);
P2graph = [P2graph, P2graphthis];
end
% y12 = P1.* (V1./x12).^y;


plot(x12,P2graph)
plot(V2,P2,'r*')
plot([V2 V3], [P2 P3])
plot(V3,P3,'r*')
plot([V3 V4], [P3 P4])
plot(V4,P4,'r*')
x45 = linspace(V4,V5,100);
P5graph = [];
for i = 1:100   
[T5graph, P5graphthis, V5graph] = JANAF(T4, P4, m, V4, 0, 0, x45(i)/V4, 4);
P5graph = [P5graph, P5graphthis];
end

% y45 = P4.* (V4./x45).^y;
plot(x45,P5graph)
plot(V5,P5,'r*')
plot([V5 V6], [P5 P6])
plot(V6,P6,'r*')
plot([V6 V7], [P6 P7])
plot(V5,P62,'r*')
plot(V7,P7,'r*')
plot([V7 V8], [P7 P8])
plot(V61,P61,'r*')



x561 = linspace(V5,V61,1000);
y561 = P5.* (V5./x561).^y;
plot(x561,y561,'--')
plot([V61 V6], [P61 P6],'--')
% patch([x561 fliplr(x561)], [y561 fliplr(linspace(P6,P61,1000))], 'r')


text(V8,P8,'8','FontSize',14)
text(V1,P1,'1','FontSize',14)
text(V2,P2,'2','FontSize',14)
text(V3,P3,'3','FontSize',14)
text(V4,P4,'4','FontSize',14)
text(V5,P5,'5','FontSize',14)
text(V6,P6,'6','FontSize',14)
text(V7,P7,'7','FontSize',14)
text(V61,P61,'6''','FontSize',14)
text(V5,P62,'6''''','FontSize',14)

set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
ylim([Patm Plim+3*10^6])
% title('PV Diagram Dual Turbocharged Miller Cycle')
title(['ALT : ' num2str(Alt) '  Boost : ' num2str(Boost) '   Pmax : ' num2str(round(P4/10^6,1))... 
    ' MPa' '   Tmax : ' num2str(round(T4)) '   Epsilon : ' num2str(Eps) ...
    '   Sigma : ' num2str(Sigma) ...
    ' Cruise Power : ' num2str(round(MaxPower)) ' Hp'])
xlabel('Volume (m^3)')
ylabel('Pressure [LogScale] (Pa)')
text(4E-3,10^6,['Efficiency :  ' num2str(round(Eff,2))])
text(4E-3,8*10^5,['bmep( MPa ) :  ' num2str(round(bmep/10^6,2))])
text(4E-3,6*10^5,['BSFC( lb/hp-h ) :  ' num2str(round(bsfc_bhp,3))])
text(4E-3,4.5*10^5,['RPM(' num2str(Power) 'Hp) : ' num2str(round(RPM))])
hold off  
end
%  '  Power : ' num2str(Power) ' Hp'
%% Table
if (Table == 1)
Params  =   {'Temp','Pressure (MPa)','Volume'};
STatm   =   [Tamb; Patm / 10^6; V8];
ST1     =   [T1; P1 / 10^6; V1];
ST2     =   [T2; P2 / 10^6; V2];
ST3     =   [T3; P3 / 10^6; V3];
ST4     =   [T4; P4 / 10^6; V4];
ST5     =   [T5; P5 / 10^6; V5];
ST6     =   [T6; P6 / 10^6; V6];
ST7     =   [T7; P7 / 10^6; V7];
ST8     =   [T8; P8 / 10^6; V7];
ST61    =   [T61; P61 / 10^6; V61];
ST62    =   [T62; P62 / 10^6; V5];
table(STatm, ST8, ST1, ST2, ST3, ST4, ST5, ST6, ST7, ST61, ST62, 'RowNames',Params)
end

%% Catches
if ((Tlim-T2)<0)
    error('T2 Exceeds Maximum Temperature')
end
if ((Plim-P4) < 0)
    error('P4 Exceeds Maximum Pressure')
end
if ((Tlim-T4)<0)
    error('T4 Exceeds Maximum Temperature')
end

if (P62 < Pe)
    error('Turbocharger Cannot Provide Necessary Compression')
end
% if (Alt > 0)
% if (RPM > 5000)
%     error('Engine RPM too High, Increase MEP')
% end
% end