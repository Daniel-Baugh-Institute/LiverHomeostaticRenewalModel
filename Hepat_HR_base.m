function [t,cellh,k] = Hepat_HR_base(feedBackInput,insultInput,startBal)
% Author: Daniel Cook
% Date: 01/01/2016
% Copyrighted under Creative Commons Share Alike

% This function runs the base model presented in "Cell network modeling to study
% dynamic robust control of the liver homeostatic renewal process" by Cook,
% Ogunnaike, and Vadigepalli
% The model takes as inputs the feedback mechanisms to include in the base
% model (see "Feedback mechanisms", below), a single insult of interest
% (see "Insults", below), and the starting balance of the SR_high and SR_low cells 
% as a fraction of total liver mass. 

% To run, use: [t,cellh,k] = Hepat_HR_base([1,1,1,0],1,[.05 .95]);
% k = parameter values, R1-R4 = Robustness per disturbance, OR = overall robustness

% Set plotting and printing (1=show results, 2=suppress results)
shouldPlot = 1; colorChoice = 'k';
shouldPrint = 1;

% Feedback mechanisms
% F(1) = implicit competition (Model A)
% F(2) = product inhibition of proliferation (Model B)
% F(3) = product inhibition of transitions (Model C)
% F(4) = alternate populations (i.e. stem cells) (Model D)

% Insults
% 0 = no insult
% 1 = increased apoptosis, 2 = decreased apoptosis
% 3 = increased proliferation, 4 = decreased proliferation
% 5 = periodic apoptosis increase
% 6 = periodic apoptosis increase w/ same width
% 7 = Periodic apoptosis increase w/ differing magnitudes
% 8 = Periodic apoptosis decrease
% 9 = Sinusoidal apoptosis rate
% 10= Unit step change in apoptosis rate

%% Section 1: Set parameter values given which feedbacks are present
% Declare global variables & input values
global F y0 insult insultFrequency
F = feedBackInput;
insult = insultInput;

% Set insult frequency
% Natural frequencies = 1.4224 & 0.3086
insultFrequency = 2; % 2 = 1/2 day on, 1/2 day off
% Use blue for high frequency (>1), red for low frequency (<1)


% Set time (in days)
tStart = 0;
tEnd = 365*2.5; % Number of years to run simulation

% Table 1: Initial hepatocyte states, cell = (SR_high,SR_low)
init_SR_high = .05; % Initial level
cellh = [init_SR_high,1-init_SR_high];
x0 = cellh;
% Set steady-state conditions for changing starting points
y0 = x0;

% Change starting concentrations if startBal parameter is used
x0(1) = startBal(1);
x0(2) = startBal(2); 

% Define model parameters
k = zeros(13,1);

% Table 2: Physiological parameters
% Rates are in doublings/day
k(1) = (1/14);  % k_prol^SR_high
k(2) = .05*1.5; % K_cap^SR_high (carrying capacity of SR_high hepatocytes)
k(3) = .95*1.5; % K_cap^SR_low (carrying capacity of SR_low hepatocytes)

% Table 3: Literature and steady-state constrained parameters
k(4) = .5*k(1); % k_prol^SR_low (Proliferation rate is half of SR_high) (Lit.)

% Additional parameters for specific feedbacks
k(8) = 1;       % k_P-env: Tissue microenvironment effect SR_high proliferation
k(9) = 0.95;    % k_T-env: Tissue microenvironment effect on cell transition
k(10) = 0.1;    % k_A
k(11) = 100;    % C1
k(12) = 8;      % C2 (This starts renewing SR_high cells at ~.025)
k(13) = 8;      % kRenew

% Set feedback conditions
if F(1) == 0
    feedBack1a = 1; feedBack1b = 1;
elseif F(1) == 1
    feedBack1a = (1-cellh(1)/k(2));
    feedBack1b = (1-cellh(2)/k(3));
end
if F(2) == 0
    feedBack2 = 1;
elseif F(2) == 1
    feedBack2 = k(8)/(cellh(1)+k(10)*cellh(2));
end
if F(3) == 0
    feedBack3 = 1;
elseif F(3) == 1
    feedBack3 = (k(9))/cellh(2);
%     hillCoeff = 1;
%     feedBack3 = 2*k(9)^hillCoeff/(k(9)^hillCoeff+cellh(2)^hillCoeff); % Alternate feedback formulation
end
if F(4) == 0
    feedBack4 = 0;
elseif F(4) == 1
    feedBack4 = k(13)*(1/(1+exp(k(11)*cellh(1)+k(12))));
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
    solverfun1c = @(x)abs(k(1)*feedBack1a*x/(cellh(1)+k(10)*cellh(2)) - k(1));
    [k(8),fval] = fminsearch(solverfun1c,1);
    feedBack2 = k(8)/(cellh(1)+k(10)*cellh(2));
    solverfun1b = @(x)abs(x*feedBack1b - k(4));
    [k(4),fval] = fminsearch(solverfun1b,1);
end
if F(2) == 1 && F(1) == 1
    solverfun1a = @(x)abs(x*feedBack1a - k(1));
    [k(1),fval] = fminsearch(solverfun1a,1);
    solverfun1c = @(x)abs(k(1)*feedBack1a*x/(cellh(1)+k(10)*cellh(2)) - k(1));
    [k(8),fval] = fminsearch(solverfun1c,1);
    feedBack2 = k(8)/(cellh(1)+k(10)*cellh(2));
    solverfun1b = @(x)abs(x*feedBack1b - k(4));
    [k(4),fval] = fminsearch(solverfun1b,1);
end
% Solve for kapHigh (using the sum of equations 1 and 2)
solverfun = @(x)abs(k(1)*cellh(1)*feedBack1a*feedBack2 + feedBack4 - cellh(1)*x + k(4)*cellh(2)*feedBack1b - cellh(2)*x);
[k(5),fval] = fminsearch(solverfun,1);
k(6) = k(5); % set kapLow = kapHigh
% Solve for kT (using equation 1)
solverfun2 = @(x)abs(k(1)*cellh(1)*feedBack1a*feedBack2 + feedBack4 - cellh(1)*k(5) - cellh(1)*x*feedBack3);
[k(7),fval2] = fminsearch(solverfun2,1);


%% Section 2: Run model
% Call ODE solver
timeStep = 0.1; % Days
[t,x] = ode45(@(t,x)odefun(t,x,k,x0), [tStart:timeStep:tEnd], x0);
cellh = x(:,1:2);

%% Section 3: Plot results
% Plot results
if shouldPlot == 1
    % Hepatocyte Populations
    figure(1); hold on; plot(t,cellh(:,1),'-','color',colorChoice,'linewidth',2);
    plot(t,cellh(:,2),'--','color',colorChoice,'linewidth',2);
    set(gca,'fontsize',18,'linewidth',2); box off
    xlabel('Time (Days)'); ylabel('Hepatocyte Populations')
    legend('SR_h_i_g_h','SR_l_o_w')

    % Hepatocyte Population Phase Plane
    figure(2); hold on; plot(cellh(:,1),cellh(:,2),'-','color',colorChoice,'linewidth',2);
    set(gca,'fontsize',18,'linewidth',2); box off
    ylabel('SR_l_o_w Hepatocytes')
    xlabel('SR_h_i_g_h Hepatocytes')
        
    % If insult is periodic, plot insult period
    if insult == 5 || insult == 8
        yTime = 0:.1:tEnd;
        yInsult = zeros(1,length(yTime));
        yInsult(-sin(insultFrequency*3.14*yTime) > 0) = 1;
        figure(3); hold on; plot(yTime,yInsult,'-','color',colorChoice,'linewidth',3)
        set(gca,'fontsize',18,'linewidth',2); box off
        ylabel('Relative apoptosis Rate'); xlabel('Time (Days)')
        xlim([0 30])
    end
    
    % If insult is periodic, plot insult period
    if insult == 6
        yTime = 0:.1:tEnd;
        yInsult = zeros(1,length(yTime));
        yInsult(-sin(insultFrequency*3.14*yTime) > sin(pi/4*(2-insultFrequency))) = 1;
        figure(3); hold on; plot(yTime,yInsult,'-','color',colorChoice,'linewidth',3)
        set(gca,'fontsize',18,'linewidth',2); box off
        ylabel('Relative apoptosis Rate'); xlabel('Time (Days)')
        xlim([0 30])
    end
    
    % If insult is periodic, plot insult period
    if insult == 7
        yTime = 0:.1:tEnd;
        yInsult = zeros(1,length(yTime));
        yInsult(-sin(insultFrequency*3.14*yTime) > 0) = 2;
        figure(3); hold on; plot(yTime,yInsult,'-','color',colorChoice,'linewidth',3)
        set(gca,'fontsize',18,'linewidth',2); box off
        ylabel('Relative apoptosis Rate'); xlabel('Time (Days)')
        xlim([0 30])
    end
    
    % If insult is sinusoidal, plot insult period
    if insult == 9 || insult == 10
        yTime = 0:.01:tEnd;
        yInsult = sin(insultFrequency*3.14*yTime)*k(6) + k(6);
        figure(3); hold on; plot(yTime,yInsult,'-','color',colorChoice,'linewidth',3)
        set(gca,'fontsize',18,'linewidth',2); box off
        ylabel('Apoptosis Rate'); xlabel('Time (Days)')
        xlim([0 30])
    end
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
SrHigh = x(1); SrLow = x(2);
SrHighIC = y0(1); SrLowIC = y0(2);
kprolHigh = k(1); kCapHigh = k(2); kapHigh = k(5);
kCapLow = k(3); kprolLow = k(4); kapLow = k(6); 
kT = k(7);
kPenv = k(8); kTenv = k(9); kA = k(10);
kRenew = k(13); C1 = k(11); C2 = k(12);

% Set insults
% Insult 1: Increased apoptosis rate (50% for 30 days)
if insult == 1
    if t >= 10 && t <= 40
        kapLow = 1.5*k(6);
        kapHigh = 1.5*k(5);
    end
end

% Insult 2: Decreased apoptosis rate (50% for 30 days)
if insult == 2
    if t >= 10 && t <= 40
        kapLow = 0.5*k(6);
        kapHigh = 0.5*k(5);
    end
end

% Insult 3: Increased proliferation rate (50% for 30 days)
if insult == 3
    if t >= 10 && t <= 40
        kprolLow = 1.5*k(4);
        kprolHigh = 1.5*k(1);
    end
end

% Insult 4: Decreased proliferation rate (50% for 30 days)
if insult == 4
    if t >= 10 && t <= 40
        kprolLow = 0.5*k(4);
        kprolHigh = 0.5*k(1);
    end
end

% Insult 5: Periodic apoptosis increase
if insult == 5
    if -sin(insultFrequency*3.14*t) > 0
        kapLow = 1.5*k(6);
        kapHigh = 1.5*k(5);
    end
end

% Insult 6: Periodic apoptosis increase w/ same width
if insult == 6
    if -sin(insultFrequency*3.14*t) > sin(pi/4*(2-insultFrequency))
        kapLow = 1.5*k(6);
        kapHigh = 1.5*k(5);
    end
end

% Insult 7: Periodic apoptosis increase w/ differing magnitudes
if insult == 7
    if -sin(insultFrequency*3.14*t) > 0
        kapLow = (0.25)*k(6)+k(6);
        kapHigh = (0.25)*k(5)+k(5);
    end
end

% Insult 8: Periodic apoptosis decrease
if insult == 8
    if -sin(insultFrequency*3.14*t) > 0
        kapLow = 0.5*k(6);
        kapHigh = 0.5*k(5);
    end
end

% Insult 9: Sinusoidal apoptosis rate
if insult == 9
    kapLow = sin(insultFrequency*3.14*t)*k(6) + k(6);
    kapHigh = sin(insultFrequency*3.14*t)*k(5) + k(6);
end

% Insult 10: Unit step change in apoptosis rate
if insult == 10
    if t >= 10
        kapLow = 0.5*k(6)+k(6);
        kapHigh = 0.5*k(5)+k(5);
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
    feedBack1a = (1-SrHigh/kCapHigh);
    feedBack1b = (1-SrLow/kCapLow);
end
if F(2) == 0
    feedBack2 = 1;
elseif F(2) == 1
    feedBack2 = kPenv/(SrHigh + kA*SrLow);
end
if F(3) == 0
    feedBack3 = 1;
elseif F(3) == 1
    feedBack3 = kTenv/SrLow;
%     hillCoeff = 1;
%     feedBack3 = 2*kTenv^hillCoeff/(kTenv^hillCoeff+A2Minus^hillCoeff); % Alternate feedback
end
if F(4) == 0
    feedBack4 = 0;
elseif F(4) == 1
    feedBack4 = kRenew*(1/(1+exp(C1*SrHigh+C2)));
end

% Equations
r1 = SrHigh*kprolHigh*feedBack1a*feedBack2 + feedBack4; % SR_high growth
r2 = SrHigh*kT*feedBack3; % Transition from SR_high to SR_low
r3 = SrHigh*kapHigh; % SR_high Hepatocyte apoptosis
r4 = SrLow*kprolLow*feedBack1b; % SR_low growth
r5 = SrLow*kapLow; % SR_low Hepatocyte apoptosis

% Differential equations
dxdt(1) = r1 - r2 - r3; % SR_high Hepatocytes
dxdt(2) = r4 + r2 - r5; % SR_low Hepatocytes

dxdt = dxdt';
end