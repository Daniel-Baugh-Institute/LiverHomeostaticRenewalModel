function [t,cellh,k] = Hepat_HR_senes(feedBackInput,insultInput)
% Author: Daniel Cook
% Date: 01/01/2016
% Copyrighted under Creative Commons Share Alike
% To run, use: [t,cellh,k] = Hepat_HR_senes([1,1,1,0],1);
% k = parameter values,

% This function runs the hepatocyte senescense model presented in "Cell network modeling to study
% dynamic robust control of the liver homeostatic renewal process" by Cook, 
% Ogunnaike, and Vadigepalli. 
% The model takes as inputs the feedback mechanisms to include in the base
% model (see "Feedback mechanisms", below) and a single insult of interest
% (see "Insults", below). The initial fractions of SR_high and SR_low cells 
% are fixed at 5% and 95% of liver mass, respectively. The percentage of 
% hepatocytes that become senescent is fixed at 99%, but can be changed in
% the code below.

% Set plotting and printing (1=show results, 2=suppress results)
shouldPlot = 1; colorChoice = 'k';
shouldPrint = 2;

% Feedback mechanisms
% F(1) = implicit competition (Model A)
% F(2) = product inhibition of proliferation (Model B)
% F(3) = product inhibition of transitions (Model C)
% F(4) = alternate populations (i.e. stem cells) (Model D)

% Insults
% 0 = no insult
% 1 = transient increased apoptosis, 
% 2 = sustained increased apoptosis apoptosis
% 3 = Removal of 40% of senenscent hepatocytes

%% Section 1: Set parameter values given which feedbacks are present
% Declare global variables & input values
global F y0 insult insultFrequency
F = feedBackInput;
insult = insultInput;

% Set percentage of hepatocytes that become senescent
senes = 0.99;

% Set time (in days)
tStart = 0;
tEnd = 365*1; % Number of years to run simulation

% Table 1: Initial hepatocyte states, cell = (SR_high,SR_low,senescent)
init_SR_high = .05; % Initial level
init_SR_low = 1-init_SR_high; % Initial level
cellh = [init_SR_high*(1-senes),init_SR_low*(1-senes),senes];
x0 = cellh;
% Set steady-state conditions for changing starting points
y0 = x0;

% Define model parameters
k = zeros(13,1);

% Insult 3: Removal of 40 percent of senescent hepatocytes
if insult == 3
    x0(3) = (1-.4)*senes;
end

% Table 2: Physiological parameters
% Rates are in doublings/day
k(1) = (1/14); % k_prol^SR_high
k(2) = .05*1.5; %0.075; % K_cap^SR_high (carrying capacity of SR_high hepatocytes)
k(3) = .95*1.5; %0.975; % K_cap^SR_low (carrying capacity of SR_low hepatocytes)

% Table 3: Literature and steady-state constrained parameters
k(4) = .5*k(1); % k_prol^SR_low (Proliferation rate is half of SR_low) (Lit.)

% Additional parameters for specific feedbacks
k(8) = 1; % k_P-env: Tissue microenvironment effect SR_high proliferation
k(9) = 0.95; % k_T-env: Tissue microenvironment effect on cell transition
k(10) = 0.1; % k_A
k(11) = 100; % C1
k(12) = 8; % C2 (This starts renewing SR_high cells at ~.025)
k(13) = 8; % kRenew

% Set feedback conditions
if F(1) == 0
    feedBack1a = 1; feedBack1b = 1;
elseif F(1) == 1
    feedBack1a = (1-init_SR_high/k(2));
    feedBack1b = (1-init_SR_low/k(3));
end
if F(2) == 0
    feedBack2 = 1;
elseif F(2) == 1
    feedBack2 = k(8)/(init_SR_high+k(10)*init_SR_low);
end
if F(3) == 0
    feedBack3 = 1;
elseif F(3) == 1
    feedBack3 = (k(9))/init_SR_low;
end
if F(4) == 0
    feedBack4 = 0;
elseif F(4) == 1
    feedBack4 = k(13)*(1/(1+exp(k(11)*init_SR_high+k(12))));
end

% Set steady-state constrained parameters
% Set effective proliferation rate equal to observed proliferation rate
if F(1) == 1 && F(2) == 0
    solverfun1a = @(x)abs(x*feedBack1a*feedBack2 - k(1));
    [k(1),fval] = fminsearch(solverfun1a,1);
    solverfun1b = @(x)abs(x*feedBack1b - k(4));
    [k(4),fval] = fminsearch(solverfun1b,1);
end
if F(2) == 1 && F(1) == 0
    solverfun1c = @(x)abs(k(1)*feedBack1a*x/(init_SR_high+k(10)*init_SR_low) - k(1));
    [k(8),fval] = fminsearch(solverfun1c,1);
    feedBack2 = k(8)/(init_SR_high+k(10)*init_SR_low);
    solverfun1b = @(x)abs(x*feedBack1b - k(4));
    [k(4),fval] = fminsearch(solverfun1b,1);
end
if F(2) == 1 && F(1) == 1
    solverfun1a = @(x)abs(x*feedBack1a - k(1));
    [k(1),fval] = fminsearch(solverfun1a,1);
    solverfun1c = @(x)abs(k(1)*feedBack1a*x/(init_SR_high+k(10)*init_SR_low) - k(1));
    [k(8),fval] = fminsearch(solverfun1c,1);
    feedBack2 = k(8)/(init_SR_high+k(10)*init_SR_low);
    solverfun1b = @(x)abs(x*feedBack1b - k(4));
    [k(4),fval] = fminsearch(solverfun1b,1);
end
% Solve for kapHigh (using the sum of equations 1 and 2)
solverfun = @(x)abs(k(1)*init_SR_high*feedBack1a*feedBack2 + feedBack4 - init_SR_high*x + k(4)*init_SR_low*feedBack1b - init_SR_low*x);
[k(5),fval] = fminsearch(solverfun,1);
k(6) = k(5); % set kapLow = kapHigh
% Solve for kT (using equation 1)
solverfun2 = @(x)abs(k(1)*init_SR_high*feedBack1a*feedBack2 + feedBack4 - init_SR_high*k(5) - init_SR_high*x*feedBack3);
[k(7),fval2] = fminsearch(solverfun2,1);

%% Section 2: Run model
% Call ODE solver
timeStep = 0.1; % Days
[t,x] = ode15s(@(t,x)odefun(t,x,k,x0), [tStart:timeStep:tEnd], x0);
cellh = x(:,1:3);

%% Section 3: Plot results
% Plot results
if shouldPlot == 1
    % Hepatocyte Populations
    figure(1); hold on; plot(t,cellh(:,1),'-','color',colorChoice,'linewidth',2);
    plot(t,cellh(:,2),'--','color',colorChoice,'linewidth',2);
    plot(t,cellh(:,3),'-.','color',colorChoice,'linewidth',2);
    set(gca,'fontsize',18,'linewidth',2); box off
    xlabel('Time (Days)'); ylabel('Hepatocyte Populations')
    legend('SR_h_i_g_h','SR_l_o_w','Senescent')

    % Hepatocyte Population Phase Plane
    figure(2); hold on; plot(cellh(:,1),cellh(:,2),'-','color',colorChoice,'linewidth',2);
    set(gca,'fontsize',18,'linewidth',2); box off
    ylabel('SR_l_o_w Hepatocytes')
    xlabel('SR_h_i_g_h Hepatocytes')
    
    timepoints = [0 10 40 60 90 120 150 180 210 240 270 300 330 365]; % 1 year
     for i=1:length(timepoints)
        plotPoint(i) = find(t >= timepoints(i),1);
     end
     figure(2); hold on; plot(cellh(plotPoint,1),cellh(plotPoint,2),'o','color',colorChoice,'markerFaceColor',colorChoice)
    
    
    % Total liver mass
    figure(3); hold on; plot(t,sum(cellh'),'-','color',colorChoice,'linewidth',2)
    set(gca,'fontsize',18,'linewidth',2); box off
    xlabel('Time (Days)')
    ylabel('Relative Liver Mass')
    
end

%% Section 4: Print results
% Calculate recovery time
% y0 is the steady-state value for cellh
recoveryTime = 1;
temp = length(find(cellh((40/timeStep+1):end,1)./y0(1)>=.99 & cellh((40/timeStep+1):end,1)./y0(1)<= 1.1,1));
if temp > 0
    temp = length(find(cellh((40/timeStep+1):end,2)./y0(2)>=.99 & cellh((40/timeStep+1):end,2)./y0(2)<= 1.1,1));
    if temp > 0 
        tRecovery(1) = t((40/timeStep+1)+find(cellh((40/timeStep+1):end,1)./y0(1)>=.99 & cellh((40/timeStep+1):end,1)./y0(1)<= 1.1,1))-40;
        tRecovery(2) = t((40/timeStep+1)+find(cellh((40/timeStep+1):end,2)./y0(2)>=.99 & cellh((40/timeStep+1):end,2)./y0(2)<= 1.1,1))-40;
        recoveryTime = max(tRecovery);
    end
end

% Calculate overall deviation
deviation = zeros(size(cellh,1),size(cellh,2));
for i=1:size(cellh,1)
    deviation(i,:) = abs(cellh(i,:) - y0)*timeStep;
end
sumDeviation = sum(deviation);

% Calculate robustness
if insult == 0
    R = 'No insult';
elseif insult ~= 0
    R = sumDeviation(1)*sumDeviation(2)*recoveryTime;
end

% Print results
if shouldPrint == 1
    fprintf('\n Robustness metric score follows: \n');
    R
end
end

function dxdt = odefun(t,x,k,x0)
% Set global parameters
global F y0 insult insultFrequency

% Set aliases
SrHigh = x(1); SrLow = x(2); Senescent = x(3);
SrHighIC = y0(1); SrLowIC = y0(2);
kprolHigh = k(1); kCapHigh = k(2); kapHigh = k(5);
kCapLow = k(3); kprolLow = k(4); kapLow = k(6); 
kT = k(7);
kPenv = k(8); kTenv = k(9); kA = k(10);
kRenew = k(13); C1 = k(11); C2 = k(12);

% Set insults
% Insult 1: Sustained increased apoptosis rate (50% for 30 days)
if insult == 1
    if t >= 10 && t <= 40
        kapLow = 1.5*k(6);
        kapHigh = 1.5*k(5);
    end
end

% Insult 2: Sustained increased apoptosis rate (50%)
if insult == 2
    if t >= 10
        kapLow = 1.5*k(6);
        kapHigh = 1.5*k(5);
    end
end

% Feedback mechanisms
% F(1) = implicit competition (Model A)
% F(2) = product inhibition of proliferation (Model B)
% F(3) = product inhibition of transitions (Model C)
% F(4) = alternate populations (i.e. stem cells) (Model D)

% Set feedback conditions
if F(1) == 0
    feedBack1a = 1; feedBack1b = 1;
elseif F(1) == 1
    feedBack1a = (1-(SrHigh+.05*Senescent)/kCapHigh);
    feedBack1b = (1-(SrLow+.95*Senescent)/kCapLow);
end
if F(2) == 0
    feedBack2 = 1;
elseif F(2) == 1
    feedBack2 = kPenv/((SrHigh+.05*Senescent) + kA*(SrLow+.95*Senescent));
end
if F(3) == 0
    feedBack3 = 1;
elseif F(3) == 1
    feedBack3 = kTenv/(SrLow+.95*Senescent);
end
if F(4) == 0
    feedBack4 = 0;
elseif F(4) == 1
    feedBack4 = kRenew*(1/(1+exp(C1*(SrHigh+.05*Senescent)+C2)));
end

% Equations
r1 = SrHigh*kprolHigh*feedBack1a*feedBack2 + feedBack4; % SrHigh growth
r2 = SrHigh*kT*feedBack3; % Transition from SrHigh to SrLow
r3 = SrHigh*kapHigh; % SrHigh Hepatocyte apoptosis
r4 = SrLow*kprolLow*feedBack1b; % SrLow growth
r5 = SrLow*kapLow; % SrLow Hepatocyte apoptosis

r6 = Senescent*kapHigh;

% Differential equations
dxdt(1) = r1 - r2 - r3; % A2+ Hepatocytes
dxdt(2) = r4 + r2 - r5; % A2- Hepatocytes
dxdt(3) = -r6; % Senescent hepatocytes

dxdt = dxdt';
end