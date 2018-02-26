%% Most Relevant Design Factors (NEED TO ENSURE NetTurbo > 0!!!)
% Diesel  =   2;          % Set to 0 if Dual, 1 if diesel, 2 if avgas

%% Input Block
Power   =   732;        % Horsepower
Alt     =   6096;       % m 6096 is 20,000
t_comp_big = 5.25;      % Time for comparing fuel consumption for bigger 
t_comp_small=6.5;       % Time for comparing fuel consumption for smaller 
% Sigma   =   1.3;          % 1-1' Additional Expansion Ratio [] 1.4
% sigma_high = 3.5;

for k = 1:2
    Diesel = k;
    if (k==1)
        sigma_high = 1.9; end
    if (k==2)
        sigma_high = 3.1; end

trials = round(10*(sigma_high - 1) + 1);
sigmas = linspace(1,sigma_high,trials);
VDS = zeros(1,trials);
BSFCS = zeros(1,trials);

BoostParams = zeros(1,trials);
VdParams = zeros(1,trials);
for j = 1:trials
Sigma = sigmas(j);
% This all gets reassigned unless Dual Cycle Selected
Boost   =   3.5;        % Ratio of Intake to Atm Pressure(Pi/Patm) 1.1254
if(Alt == 0)Boost = 0.66*Boost;end % Use WASTEGATE to only have boost at alt
if(Boost < 1) Boost = 1; end; % just in case ratio drops it too low
xq      =   1;        % Fraction of Heat Input that is ISOCHORIC
Eps     =   22;          % Effective Compression Ratio [] i like 20
Vd      =   1;          % Displacement volume (Liters)
Lean    =   0.5;       % Fraction of stoich fuel to air ratio, .5 gud, 0.606 max for diesel
%% NOTE!! THESE ARE ALSO HARDCODED INTO THE TRIAL CODE!!!!
Plim = 1.034E7; 
Tlim    =   4000; 

if (Diesel == 1) 
    xq = 0.2; 
    Plim = 2E7;
    P2lim = 2E7;
    Lean = 0.606;
    Eps = 20;
end
if (Diesel == 2)
    xq = 0.8;
    Plim = 1.034E7;
    P2lim = 0.246E7;
    Lean = 0.89;
    Eps = 10;
end


%% Boost As High as Possible, Don't Modify

Boosts = linspace(1,5,41);
for i = 1:40
    
    [T2, T4, P4, P62, Pe, MaxPower, Eff, bsfc_bhp] = ...
        MillerV3JANAF_trial(0, Diesel, Power, Alt, Boosts(i), xq, Eps, ...
        Sigma, Vd, Lean);
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % TESTING FOR COMPRESSION IN AVGAS?
    if (Diesel == 2)
    if ((753-T2) < 0) % its 748 @ T2 for C ratio 12 (see sources) at no boost, no sigma, alt=0, so good appx limit
%     disp('T2 Temp Limit')
%     close all; 
    
    [T2, T4, P4, P62, Pe, MaxPower, Eff, bsfc_bhp] = ...
        MillerV3JANAF_trial(0, Diesel, Power, Alt, Boosts(i-1), xq, Eps, ...
        Sigma, Vd, Lean);
    break
    end
    end
    %END AVGAS TESTING
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    if ((Tlim-T2)<0)
%     disp('T2 Exceeds Maximum Temperature')
%     close all; 
    [T2, T4, P4, P62, Pe, MaxPower, Eff, bsfc_bhp] = ...
        MillerV3JANAF_trial(0, Diesel, Power, Alt, Boosts(i-1), xq, Eps, ...
        Sigma, Vd, Lean); 
    break
    end
    if ((Plim-P4) < 0)
%     disp('P4 Exceeds Maximum Pressure')
%     close all; 
    
    [T2, T4, P4, P62, Pe, MaxPower, Eff, bsfc_bhp] = ...
        MillerV3JANAF_trial(0, Diesel, Power, Alt, Boosts(i-1), xq, Eps, ...
        Sigma, Vd, Lean);
    break
    end
    if ((Tlim-T4)<0)
%     disp('T4 Exceeds Maximum Temperature')
%     close all; 
    [T2, T4, P4, P62, Pe, MaxPower, Eff, bsfc_bhp] = ...
        MillerV3JANAF_trial(0, Diesel, Power, Alt, Boosts(i-1), xq, Eps, ...
        Sigma, Vd, Lean);
    break
    end
    if (P62 < Pe)
%     disp('Turbocharger Cannot Provide Necessary Compression')
%     close all; 
    [T2, T4, P4, P62, Pe, MaxPower, Eff, bsfc_bhp] = ...
        MillerV3JANAF_trial(0, Diesel, Power, Alt, Boosts(i-1), xq, Eps, ...
        Sigma, Vd, Lean);
    break
    end
end
Boost = Boosts(i-1);
% disp(Eff)
% disp(MaxPower)
Vd = (Power/MaxPower)*Vd;
Vd = ceil(Vd*10)/10;
% while (MaxPower < Power)
%     Vd = Vd+0.1;
%     [T2, T4, P4, P62, Pe, MaxPower, Eff, bsfc_bhp] = ...
%         MillerV3JANAF_trial(0, Diesel, Power, Alt, Boost, xq, Eps, ...
%         Sigma, Vd, Lean);
% end
VDS(j) = Vd;
BSFCS(j) =  bsfc_bhp;
type = strcat(num2str(j), ' of ',' ', num2str(trials), ' trials completed.');
disp(type)
% Eff
% bsfc_bhp
figure(k)
[T2, T4, P4, P62, Pe, MaxPower, Eff, bsfc_bhp] = ...
        MillerV3JANAF_trial(1, Diesel, Power, Alt, Boost, xq, Eps, ...
        Sigma, Vd, Lean); 
BoostParams(j) = Boost;   
end
params = [sigmas', VDS', BSFCS', BoostParams'];
if (Diesel == 1)
    DPARAMS = params;
end
if (Diesel == 2)
    AVPARAMS = params;
end
figure(3)
hold on
plot(VDS, BSFCS)
figure(4)
hold on
plot(BSFCS, VDS)
% figure(5)
% plot(sigmas, VDS)
% figure(6)
% plot(sigmas, BSFCS)
end
%%
DLENGTH = length(DPARAMS(:,2));
DBIG =  zeros(1,DLENGTH);
DSMALL = zeros(1,DLENGTH);
DVDS = DPARAMS(:,2);
DBSFCS = DPARAMS(:,3);
for q = 1:DLENGTH
    DBIG(q) = 50*DVDS(q) + DBSFCS(q)*Power*t_comp_big*2; % Times 2 means emissions get weight of fuel
    DSMALL(q) = 50*DVDS(q) + DBSFCS(q)*Power*t_comp_small*2;
end
AVLENGTH = length(AVPARAMS(:,2));
AVBIG =  zeros(1,AVLENGTH);
AVSMALL = zeros(1,AVLENGTH);
AVVDS = AVPARAMS(:,2);
AVBSFCS = AVPARAMS(:,3);
for q = 1:AVLENGTH
    AVBIG(q) = 50*AVVDS(q) + AVBSFCS(q)*Power*t_comp_big*2; % Times 2 means emissions get weight of fuel
    AVSMALL(q) = 50*AVVDS(q) + AVBSFCS(q)*Power*t_comp_small*2;
end
figure
hold on
plot(DPARAMS(:,2), DBIG)
plot(AVPARAMS(:,2), AVBIG)
plot(DPARAMS(:,2), DSMALL)
plot(AVPARAMS(:,2), AVSMALL)
xlabel('Volume (L)')
ylabel('Weight of Engine Plus Fuel')
title(['Power: ' num2str(Power) 'Aircraft Weight Cost Trade Volume'])
legend('Diesel Big','AVGas Big', 'Diesel Small', 'Avgas Small')
hold off
figure
hold on
plot(DPARAMS(:,3), DBIG)
plot(AVPARAMS(:,3), AVBIG)
plot(DPARAMS(:,3), DSMALL)
plot(AVPARAMS(:,3), AVSMALL)
xlabel('BSFC')
ylabel('Weight of Engine Plus Fuel')
title(['Power: ' num2str(Power) 'Aircraft Weight Cost Trade BSFC'])
legend('Diesel Big','AVGas Big', 'Diesel Small', 'Avgas Small')
hold off
figure
hold on
plot(DPARAMS(:,1), DBIG)
plot(AVPARAMS(:,1), AVBIG)
plot(DPARAMS(:,1), DSMALL)
plot(AVPARAMS(:,1), AVSMALL)
xlabel('SIGMA')
ylabel('Weight of Engine Plus Fuel')
title(['Power: ' num2str(Power) 'Aircraft Weight Cost Trade SIGMA'])
legend('Diesel Big','AVGas Big', 'Diesel Small', 'Avgas Small')
hold off
figure
hold on
plot(DPARAMS(:,4), DBIG)
plot(AVPARAMS(:,4), AVBIG)
plot(DPARAMS(:,4), DSMALL)
plot(AVPARAMS(:,4), AVSMALL)
xlabel('Boost')
ylabel('Weight of Engine Plus Fuel')
title(['Power: ' num2str(Power) 'Aircraft Weight Cost Trade Boost'])
legend('Diesel Big','AVGas Big', 'Diesel Small', 'Avgas Small')
hold off




