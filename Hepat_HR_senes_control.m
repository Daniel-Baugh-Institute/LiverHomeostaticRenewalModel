function [t,cellh,controlAction,k,dxdt] = Hepat_HR_senes_control(feedBackInput,insultInput,senesPct)
% Author: Daniel Cook
% Date: 01/01/2016
% Copyrighted under Creative Commons Share Alike
% To run, use: [t,cellh,controlAction,k] = Hepat_HR_senes_control([1,1,1,0],1,0.99);
% k = parameter values, R1-R4 = Robustness per disturbance, OR = overall robustness

% This function runs the hepatocyte senescense model with non-parenchymal cell
% control presented in "Elucidating the mechanisms of 
% dynamic and robust control of the liver homeostatic renewal process: Cell network 
% modeling and analysis" by Cook, Manchel, Ogunnaike, and Vadigepalli.

% The model takes as inputs the feedback mechanisms to include in the base
% model (see "Feedback mechanisms", below), a single insult of interest
% (see "Insults", below), and the starting balance of the SR_high and SR_low cells 
% as a fraction of total liver mass. This function includes control of
% hepatocyte levels by non-parenchymal cells. For the model without
% control, see Supplemental File S2 (Hepat_HR_base.m).

% Set plotting and printing (1=show results, 2=suppress results)
shouldPlot = 0; colorChoice = 'k';
shouldPrint = 2;

% Feedback mechanisms
% F(1) = implicit competition (Model A)
% F(2) = product inhibition of proliferation (Model B)
% F(3) = product inhibition of transitions (Model C)
% F(4) = alternate populations (i.e. stem cells) (Model D)

% Insults
% 0 = no insult
% 1 = increased apoptosis, 2 = decreased apoptosis
% 3 = increased proliferation, 4 = decreased proliferation
% 5 = periodic apoptosis increase, 6 = periodic apoptosis decrease
% 7 = periodic proliferation increase, 8 = periodic proliferation increase
% 9 = sinusoidal apoptosis profile, 10 = sinusoidal proliferation profile
% 11 = Removal of 40% of senenscent hepatocytes

%% Section 1: Set parameter values given which feedbacks are present
% Declare global variables & input values
global F y0 insult insultFrequency
F = feedBackInput;
insult = insultInput;

% Set percentage of hepatocytes that become senescent
senes = senesPct;

% Set insult frequency
% Natural frequencies = 1.4224 & 0.3086
insultFrequency = 1.0; % 2 = 1/2 day on, 1/2 day off
% Use blue for high frequency (>1), red for low frequency (<1)


% Set time (in days)
tStart = 0;
tEnd = 365*1; % Number of years to run simulation

% Table 1: Initial hepatocyte states, cell = (SR-high,SR-low)
init_A2plus = .05; % Initial level
init_A2minus = 0.95;
cellh = [init_A2plus*(1-senes),init_A2minus*(1-senes),senes];
x0 = cellh;
% Set steady-state conditions for initialization of controller set points
y0 = [init_A2plus,init_A2minus,0];

% Set control action to 0 at baseline
x0(4:5) = [0 0]; 

% Define model parameters
k = zeros(13,1);

    % Insult 11: Removal of 40 percent of senescent hepatocytes
    if insult == 11
        x0(3) = (1-.4)*senes;
    end


% Table 2: Physiological parameters
% Rates are in doublings/day
k(1) = (1/14); % k_prol^SR-high
k(2) = .05*1.5; %0.075; % K_cap^SR-high (carrying capacity of SR-high hepatocytes)
k(3) = .95*1.5; %0.975; % K_cap^SR-low (carrying capacity of SR-low hepatocytes)

% Table 3: Literature and steady-state constrained parameters
k(4) = .5*k(1); % k_prol^SR-low (Proliferation rate is half of SR-low) (Lit.)

% Additional parameters for specific feedbacks
k(8) = 1; % k_P-env: Tissue microenvironment effect SR-high proliferation
k(9) = 0.95; % k_T-env: Tissue microenvironment effect on cell transition
k(10) = 0.1; % k_A
k(11) = 100; % C1
k(12) = 8; % C2 (This starts renewing SR-high cells at ~.025)
k(13) = 8; % kRenew

% Set feedback conditions
if F(1) == 0
    feedBack1a = 1; feedBack1b = 1;
elseif F(1) == 1
    feedBack1a = (1-init_A2plus/k(2));
    feedBack1b = (1-init_A2minus/k(3));
end
if F(2) == 0
    feedBack2 = 1;
elseif F(2) == 1
    feedBack2 = k(8)/(init_A2plus+k(10)*init_A2minus);
end
if F(3) == 0
    feedBack3 = 1;
elseif F(3) == 1
    feedBack3 = (k(9))/init_A2minus;
%     hillCoeff = 1;
%     feedBack3 = 2*k(9)^hillCoeff/(k(9)^hillCoeff+init_A2minus^hillCoeff); % Alternate feedback
end
if F(4) == 0
    feedBack4 = 0;
elseif F(4) == 1
    feedBack4 = k(13)*(1/(1+exp(k(11)*init_A2plus+k(12))));
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
    solverfun1c = @(x)abs(k(1)*feedBack1a*x/(init_A2plus+k(10)*init_A2minus) - k(1));
    [k(8),fval] = fminsearch(solverfun1c,1);
    feedBack2 = k(8)/(init_A2plus+k(10)*init_A2minus);
    solverfun1b = @(x)abs(x*feedBack1b - k(4));
    [k(4),fval] = fminsearch(solverfun1b,1);
end
if F(2) == 1 && F(1) == 1
    solverfun1a = @(x)abs(x*feedBack1a - k(1));
    [k(1),fval] = fminsearch(solverfun1a,1);
    solverfun1c = @(x)abs(k(1)*feedBack1a*x/(init_A2plus+k(10)*init_A2minus) - k(1));
    [k(8),fval] = fminsearch(solverfun1c,1);
    feedBack2 = k(8)/(init_A2plus+k(10)*init_A2minus);
    solverfun1b = @(x)abs(x*feedBack1b - k(4));
    [k(4),fval] = fminsearch(solverfun1b,1);
end
% Solve for kapHigh (using the sum of equations 1 and 2)
solverfun = @(x)abs(k(1)*init_A2plus*feedBack1a*feedBack2 + feedBack4 - init_A2plus*x + k(4)*init_A2minus*feedBack1b - init_A2minus*x);
[k(5),fval] = fminsearch(solverfun,1);
k(6) = k(5); % set kapLow = kapHigh
% Solve for kT (using equation 1)
solverfun2 = @(x)abs(k(1)*init_A2plus*feedBack1a*feedBack2 + feedBack4 - init_A2plus*k(5) - init_A2plus*x*feedBack3);
[k(7),fval2] = fminsearch(solverfun2,1);

% Change any parameter values for sensitivity plot (Figure 4)
% k(2) = .05*1.1;
% k(3) = 0.95*1.5;
% k(6) = 0.5*k(6);
% k = k+k.*randn(13,1)/100;

%% Section 2: Run model
% Call ODE solver
timeStep = 0.1; % Days
[t,x] = ode15s(@(t,x)odefun(t,x,k,x0), [tStart:timeStep:tEnd], x0);
cellh = x(:,1:3);
controlAction = x(:,4:5);

%% Section 3: Plot results
% Plot results
if shouldPlot == 1
    % Hepatocyte Populations
    figure(1); hold on; plot(t,cellh(:,1),'-','color',colorChoice,'linewidth',2);
    plot(t,cellh(:,2),'--','color',colorChoice,'linewidth',2);
    plot(t,cellh(:,3),'-.','color',colorChoice,'linewidth',2);
    set(gca,'fontsize',18,'linewidth',2); box off
    xlabel('Time (Days)'); ylabel('Hepatocyte Populations')
    legend('SR_h_i_g_h','SR_l_o_w', 'Senescent')

    % Hepatocyte Population Phase Plane
    figure(2); hold on; plot(cellh(:,1),cellh(:,2),'-','color',colorChoice,'linewidth',2);
    set(gca,'fontsize',18,'linewidth',2); box off
    ylabel('SR_l_o_w Hepatocytes')
    xlabel('SR_h_i_g_h Hepatocytes')
    
%     timepoints = [0 10 40 60 90 120 150 180 210 240 270 300 330 365]; % 1 year
%     timepoints = [0 10 40 60 90 120 150 180 210 240 270 300 330 365 730 1095]; % 3 year
%     timepoints = [0 10 40 60 90 120 150 180 210 240 270 300 330 365 730 1095 ...
%         1460 1825 2190 2555 2920 3285 3650 4015 4380 4745 5110 5475]; % 15 year
%     for i=1:length(timepoints)
%         plotPoint(i) = find(t >= timepoints(i),1);
%     end
%     figure(2); hold on; plot(cellh(plotPoint,1),cellh(plotPoint,2),'o','color',colorChoice,'markerFaceColor',colorChoice)
    
    % If insult is periodic, plot insult period
    if insult == 5 || insult == 8
        yTime = 0:.1:tEnd;
        yInsult = zeros(1,length(yTime));
        yInsult(-sin(insultFrequency*3.14*yTime) > 0) = 1;
        figure(3); hold on; plot(yTime,yInsult,'-','color',colorChoice,'linewidth',3)
        set(gca,'fontsize',18,'linewidth',2); box off
        ylabel('Relative cell death rate'); xlabel('Time (Days)')
        xlim([0 30])
    end
    
    % If insult is periodic, plot insult period
    if insult == 6
        yTime = 0:.1:tEnd;
        yInsult = zeros(1,length(yTime));
        yInsult(-sin(insultFrequency*3.14*yTime) > sin(pi/4*(2-insultFrequency))) = 1;
        figure(3); hold on; plot(yTime,yInsult,'-','color',colorChoice,'linewidth',3)
        set(gca,'fontsize',18,'linewidth',2); box off
        ylabel('Relative cell death rate'); xlabel('Time (Days)')
        xlim([0 30])
    end
    
    % If insult is periodic, plot insult period
    if insult == 7
        yTime = 0:.1:tEnd;
        yInsult = zeros(1,length(yTime));
        yInsult(-sin(insultFrequency*3.14*yTime) > 0) = 2;
        figure(3); hold on; plot(yTime,yInsult,'-','color',colorChoice,'linewidth',3)
        set(gca,'fontsize',18,'linewidth',2); box off
        ylabel('Relative cell death rate'); xlabel('Time (Days)')
        xlim([0 30])
    end
    
    % If insult is sinusoidal, plot insult period
    if insult == 9
        yTime = 0:.01:tEnd;
        yInsult = sin(insultFrequency*3.14*yTime)*k(6) + k(6);
        figure(3); hold on; plot(yTime,yInsult,'-','color',colorChoice,'linewidth',3)
        set(gca,'fontsize',18,'linewidth',2); box off
        ylabel('Cell death rate'); xlabel('Time (Days)')
        xlim([0 30])
    end


%     % Plot combined figure
    yTime = 0:.1:tEnd;
    yInsult = zeros(1,length(yTime));
    yInsult(-sin(insultFrequency*3.14*yTime) > 0)= 1.0; %sin(pi/4*(2-insultFrequency))) = 1;
%     yInsult(yTime >= 10) = 1; yInsult(yTime >= 40) = 0;
    figure(4); hold on;
    subplot(3,1,1); plot(yTime,yInsult,'k-','linewidth',2)
    set(gca,'fontsize',18,'linewidth',2); box off
    ylim([0 2])
    ylabel('Cell death rate'); %xlim([290 300]);
    subplot(3,1,2); hold on; plot(t,cellh(:,1),'-','color',colorChoice,'linewidth',2);
    set(gca,'fontsize',18,'linewidth',2); box off
    ylabel('SR_h_i_g_h'); %xlim([290 300]);
    subplot(3,1,3); hold on; plot(t,cellh(:,2),'-','color',colorChoice,'linewidth',2);
    set(gca,'fontsize',18,'linewidth',2); box off
    ylabel('SR_l_o_w'); %xlim([290 300]);
    xlabel('Time (Days)')
    
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
A2Plus = x(1); A2Minus = x(2); Senescent = x(3);
A2PlusIC = y0(1); A2MinusIC = y0(2);
kprolPlus = k(1); kCapPlus = k(2); kapPlus = k(5);
kCapMinus = k(3); kprolMinus = k(4); kapMinus = k(6); 
kT = k(7);
kPenv = k(8); kTenv = k(9); kA = k(10);
kRenew = k(13); C1 = k(11); C2 = k(12);

% Set controller parametersn
Kp1 = 37.5; Ki1 = 20.8;
Kp2 = 15.9; Ki2 = 2.5;

% Set insults
% Insult 1: Increased apoptosis rate (50% for 30 days)
if insult == 1
    if t >= 10 && t <= 40
        kapMinus = 1.5*k(6);
        kapPlus = 1.5*k(5);
%         kapMinus = 1.3*k(6); % For Kcap = 1.9x
%         kapPlus = 1.3*k(5); % For Kcap = 1.9x
%         kapMinus = 4.0*k(6); % For Kcap = 1.1x
%         kapPlus = 1.0*k(5); % For Kcap = 1.1x
        
    end
end

% Insult 2: Decreased apoptosis rate (50% for 30 days)
if insult == 2
    if t >= 10 && t <= 40
        kapMinus = 0.5*k(6);
        kapPlus = 0.5*k(5);
    end
end

% Insult 3: Increased proliferation rate (50% for 30 days)
if insult == 3
    if t >= 10 && t <= 40
        kprolMinus = 1.5*k(4);
        kprolPlus = 1.5*k(1);
    end
end

% Insult 4: Decreased proliferation rate (50% for 30 days)
if insult == 4
    if t >= 10 && t <= 40
        kprolMinus = 0.5*k(4);
        kprolPlus = 0.5*k(1);
    end
end

% Insult 5: Periodic apoptosis increase
if insult == 5
    if -sin(insultFrequency*3.14*t) > 0
        kapMinus = 1.5*k(6);
        kapPlus = 1.5*k(5);
    end
end

% Insult 6: Periodic apoptosis increase w/ same width
if insult == 6
    if -sin(insultFrequency*3.14*t) > sin(pi/4*(2-insultFrequency))
        kapMinus = 1.5*k(6);
        kapPlus = 1.5*k(5);
    end
end

% Insult 7: Periodic apoptosis increase w/ differing magnitudes
if insult == 7
    if -sin(insultFrequency*3.14*t) > 0
        kapMinus = 0.25*(0.5)*k(6)+k(6);
        kapPlus = 0.25*(0.5)*k(5)+k(5);
    end
end

% Insult 9: Sinusoidal apoptosis rate
if insult == 9
    kapMinus = sin(insultFrequency*3.14*t)*k(6) + k(6);
    kapPlus = sin(insultFrequency*3.14*t)*k(5) + k(6);
end

% Insult 10: Unit step change in apoptosis rate
if insult == 10
    if t >= 10
        kapMinus = 0.5*k(6)+k(6);
        kapPlus = 0.5*k(5)+k(5);
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
    feedBack1a = (1-(A2Plus+0.05*Senescent)/kCapPlus);
    feedBack1b = (1-(A2Minus+0.95*Senescent)/kCapMinus);
end
if F(2) == 0
    feedBack2 = 1;
elseif F(2) == 1
    feedBack2 = kPenv*(1+x(4))/((A2Plus+0.05*Senescent) + kA*(A2Minus+0.95*Senescent)); % change to 0.1*x(4) for 1% of CA1
end
if F(3) == 0
    feedBack3 = 1;
elseif F(3) == 1
    feedBack3 = kTenv*(1+x(5))/(A2Minus+0.95*Senescent); % change to 0.1*x(5) for 1% of CA2
%     hillCoeff = 1;
%     feedBack3 = 2*kTenv^hillCoeff/(kTenv^hillCoeff+A2Minus^hillCoeff); % Alternate feedback
end
if F(4) == 0
    feedBack4 = 0;
elseif F(4) == 1
    feedBack4 = kRenew*(1/(1+exp(C1*(A2Plus+.05*Senescent)+C2)));
end

% Equations
r1 = A2Plus*kprolPlus*feedBack1a*feedBack2 + feedBack4; % SR-high growth
r2 = A2Plus*kT*feedBack3; % Transition from SR-high to SR-low
r3 = A2Plus*kapPlus; % SR-high Hepatocyte cell death
r4 = A2Minus*kprolMinus*feedBack1b; % SR-low growth
r5 = A2Minus*kapMinus; % SR-low Hepatocyte apoptosis
r6 = Senescent*kapPlus;


% Differential equations
dxdt(1) = r1 - r2 - r3; % SR-high Hepatocytes
dxdt(2) = r4 + r2 - r5; % SR-low Hepatocytes
dxdt(3) = -r6; % Senescent hepatocytes

dxdt(4) = Kp1*(-dxdt(1)) + Ki1*(A2PlusIC-A2Plus); % Controller action (k_P_env)
dxdt(5) = Kp2*(-dxdt(2)) + Ki2*(A2MinusIC-A2Minus); % Controller action (k_T_env)


dxdt = dxdt';
end