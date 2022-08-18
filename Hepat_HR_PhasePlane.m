function [k0,k1,k2,k3,k4,k5] = Hepat_HR_PhasePlane
% Author: Daniel Cook
% Date: 01/01/2016
% Copyrighted under Creative Commons Share Alike
% To run, use: [k1,k2,k3,k4,k5] = Hepat_HR_PhasePlane;

% Define feedback mechanism
% 1 = none, 2 = implicit competition, 3 = metabolic demand
% 4 = WNT limited proliferation, 5 = alternate populations

% Plot results (1 = yes, 2 = no)
shouldPlot = 1;

%% Table 1: Initial hepatocyte states, cell = (A2+,A2-)
init_A2plus = .05; % Initial level
cellh = [init_A2plus,1-init_A2plus];
x0 = cellh;

% Define model parameters
k0 = zeros(13,1);

% Table 2: Physiological parameters
% Rates are in doublings/day
k0(1) = (1/14); % k_prol^A2+
k0(2) = 0.075; % K_cap^A2+ (carrying capacity of A2+ hepatocytes)
k0(3) = 1; % k_P (proportion of A2+ cells likely to transform)
k0(4) = 0.975; % K_cap^A2- (carrying capacity of A2- hepatocytes)

% Table 3: Literature and steady-state constrained parameters
k0(5) = .5*k0(1); % k_prol^A2- (Proliferation rate is half of A2+) (Lit.)

% Set values of parameter sets common for all feedback mechanisms
k1 = k0; k2 = k0; k3 = k0; k4 = k0; k5 = k0;

% No feedback (Condition 0)
k0(6) = (cellh(1)*k0(1) + cellh(2)*k0(5))/(cellh(1)+cellh(2)); % k_ap^A2+
k0(7) = k0(6); % k_ap^A2- (apoptosis rates are equal) (Lit.)
k0(8) = k0(1) - k0(6); % k_T

% Implicit hepatocyte competition for nutrients (Condition 1)
k1(6) = (cellh(1)*k1(1)*(1-cellh(1)/k1(2)) + cellh(2)*k1(5)*(1-cellh(2)/k1(4)))/(cellh(1)+cellh(2)); % k_ap^A2+ (Steady-state)
k1(7) = k1(6); % k_ap^A2- (apoptosis rates are equal) (Lit.)
k1(8) = (k1(1)*(1-cellh(1)/k1(2)) - k1(6))/k1(3); % k_T (Steady-state)

% Metabolic population requirement for Axin2+ cells (Condition 2)
k2(6) = (cellh(1)*k2(1) + cellh(2)*k2(5))/(cellh(1)+cellh(2)); % k_ap^A2+
k2(7) = k2(6); % k_ap^A2- (apoptosis rates are equal) (Lit.)
k2(8) = k2(1) - k2(6); % k_T (Steady-state)

% Proliferation of Axin2+ cells is dependent on WNT concentration (Condition 3)
k3(9) = 1; % Available WNT concentration
k3(10) = 0.1; % Fraction of Axin2- cells contributing to WNT diffusion
k3(6) = (cellh(1)*k3(1)*k3(9)/(cellh(1)+k3(10)*cellh(2))+ cellh(2)*k3(5))/(cellh(1)+cellh(2)); % k_ap^A2+
k3(7) = k3(6); % k_ap^A2- (apoptosis rates are equal) (Lit.)
k3(8) = k3(1)*k3(9)/(cellh(1)+k3(10)*cellh(2)) - k3(6); % k_T (Steady-state)

% Transition of Axin2+ cells is dependent on WNT concentration (Condition 4)
k4(9) = 1; % Available WNT concentration
k4(10) = 0.1; % Fraction of Axin2- cells contributing to WNT diffusion
k4(6) = (cellh(1)*k4(1) + cellh(2)*k4(5))/(cellh(1)+cellh(2)); % k_ap^A2+
k4(7) = k4(6); % k_ap^A2- (apoptosis rates are equal) (Lit.)
k4(8) = (k4(1) - k4(6))*k4(9)/(cellh(1)+k4(10)*cellh(2)); % k_T (Steady-state)

% Alternate population renewing Axin2+ cells (ex. stem cells) (Condition 5)
k5(11) = 100; % C1
k5(12) = 8; % C2
k5(13) = 8; % kRenew
k5(6) = (cellh(1)*k5(1)+ cellh(2)*k5(5) + k5(13)*(1/(1+exp(k5(11)*cellh(1)+k5(12)))))/(cellh(1)+cellh(2)); % k_ap^A2+
k5(7) = k5(6); % k_ap^A2- (apoptosis rates are equal) (Lit.)
k5(8) = k5(1) - k5(6) + k5(13)*(1/(1+exp(k5(11)*cellh(1)+k5(12))))/cellh(1); % k_T (Steady-state)

%% Find and plot nulclines
%% No feedback (Condition 0)
% Set aliases
A2PlusIC = x0(1); A2MinusIC = x0(2);
kprolPlus = k0(1); kCapPlus = k0(2); kapPlus = k0(6);
kCapMinus = k0(4); kprolMinus = k0(5); kapMinus = k0(7); 
kP = k0(3); kT = k0(8);
WNT = k0(9); kDiff = k0(10);
kRenew = k0(13); C1 = k0(11); C2 = k0(12);

% Equations
HPlus = 0:.01:2;
HMinus = -HPlus*kT/(kprolMinus - kapMinus);

% Plot results
if shouldPlot == 1
    % Plot phase plane results
    figure(); hold on;
    plot(HPlus,HMinus,'k-','linewidth',2);
    plot(zeros(length(HPlus),1),HPlus,'k-','linewidth',2);
    xlim([0 2]); ylim([0 2]);
    set(gca,'fontsize',18,'linewidth',2); box off
    ylabel('Axin2- Hepatocytes'); xlabel('Axin2+ Hepatocytes')
    title('No Feedback','fontSize',24);
    
    % Plot parameter values
    figure(); hold on; bar(k0,'FaceColor',[.5 .5 .5],'EdgeColor',[0 0 0],'linewidth',2)
    set(gca,'fontsize',18,'linewidth',2); box off; ylim([0 1.2])
    xlabel('Parameter Number'); ylabel('Parameter Value')
    title('No Feedback','fontsize',24)
end

%% Implicit hepatocyte competition for nutrients (Condition 1)    
% Set aliases
A2PlusIC = x0(1); A2MinusIC = x0(2);
kprolPlus = k1(1); kCapPlus = k1(2); kapPlus = k1(6);
kCapMinus = k1(4); kprolMinus = k1(5); kapMinus = k1(7); 
kP = k1(3); kT = k1(8);
WNT = k1(9); kDiff = k1(10);
kRenew = k1(13); C1 = k1(11); C2 = k1(12);

% Equations
HMinus = 0:.01:2;
A2Plus = (kprolPlus - kT*kP - kapPlus)*(kCapPlus/kprolPlus);

HPlus = (1/(kT*kP))*(HMinus*kapMinus - HMinus*kprolMinus.*(1-HMinus/kCapMinus));

% Plot results
if shouldPlot == 1
    % Plot phase plane results
    figure(); hold on;
    plot(A2Plus*ones(length(HMinus),1),HMinus,'k-','linewidth',2);
    plot(zeros(length(HMinus),1),HMinus,'k-','linewidth',2);
    plot(HPlus,HMinus,'k-','linewidth',2)
    xlim([-.1 .2]); ylim([0 1.2]);
    set(gca,'fontsize',18,'linewidth',2); box off
    ylabel('Axin2- Hepatocytes'); xlabel('Axin2+ Hepatocytes')
    title('Implicit Hepatocyte Competition','fontSize',24)
    
    % Plot parameter values vs. no feedback
    figure(); hold on; plot(k1,k0,'ko','markerFaceColor','k','markerSize',4)
    plot([0 1.5],[0 1.5],'k-','linewidth',1); xlim([0 1])
    set(gca,'fontsize',18,'linewidth',2); box off
    xlabel('Implicit competition parameter values')
    ylabel('No feedback parameter values')
    title('Implicit Hepatocyte Competition','fontsize',24)
end

%% Metabolic population requirement for Axin2+ cells (Condition 2)
% Set aliases
A2PlusIC = x0(1); A2MinusIC = x0(2);
kprolPlus = k2(1); kCapPlus = k2(2); kapPlus = k2(6);
kCapMinus = k2(4); kprolMinus = k2(5); kapMinus = k2(7); 
kP = k2(3); kT = k2(8);
WNT = k2(9); kDiff = k2(10);
kRenew = k2(13); C1 = k2(11); C2 = k2(12);

% Equations
HPlus = 0:.01:2;
A2Plus = A2PlusIC/(kT*kP)*(kprolPlus - kapPlus);

HMinus = -1/(kprolMinus - kapMinus)*(kT*kP*HPlus.^2)/A2PlusIC;

if shouldPlot == 1
    % Plot phase plane results
    figure(); hold on;
    plot(A2Plus*ones(length(HPlus),1),HPlus,'k-','linewidth',2);
    plot(HPlus,HMinus,'k-','linewidth',2)
    xlim([0 .1]); ylim([0 2]);
    set(gca,'fontsize',18,'linewidth',2); box off
    ylabel('Axin2- Hepatocytes'); xlabel('Axin2+ Hepatocytes')
    title('Metabolic Feedback','fontSize',24)
    
%     % Add arrows to plot
%     for A2Plus = 0:.01:.1
%         for A2Minus = 0:.1:2
%                 r1 = A2Plus*kprolPlus; % A2Plus Hepatocyte proliferation
%                 r2 = A2Plus*kT*kP*(A2Plus/A2PlusIC); % Transition from A2Plus to A2Minus
%                 r3 = A2Plus*kapPlus; % A2Plus Hepatocyte apoptosis
%                 r4 = A2Minus*kprolMinus; % A2Minus Hepatocyte proliferation
%                 r5 = A2Minus*kapMinus; % A2Minus Hepatocyte apoptosis
%                 dxdt(1) = r1 - r2 - r3; % A2+ Hepatocytes
%                 dxdt(2) = r4 + r2 - r5; % A2- Hepatocytes
%                 arrow([A2Plus,A2Minus],[A2Plus+dxdt(1),A2Minus+dxdt(2)])
%         end
%     end
    
    % Plot parameter values vs. no feedback
    figure(); hold on; plot(k2,k0,'ko','markerFaceColor','k','markerSize',4)
    plot([0 1.5],[0 1.5],'k-','linewidth',1); xlim([0 1])
    set(gca,'fontsize',18,'linewidth',2); box off
    xlabel('Implicit competition parameter values')
    ylabel('No feedback parameter values')
    title('Metabolic Feedback','fontsize',24)

end
%% Proliferation of Axin2+ cells is dependent on WNT concentration (Condition 3)
% Set aliases
A2PlusIC = x0(1); A2MinusIC = x0(2);
kprolPlus = k3(1); kCapPlus = k3(2); kapPlus = k3(6);
kCapMinus = k3(4); kprolMinus = k3(5); kapMinus = k3(7); 
kP = k3(3); kT = k3(8);
WNT = k3(9); kDiff = k3(10);
kRenew = k3(13); C1 = k3(11); C2 = k3(12);

% Equations
y = 0:.1:2; % Dummy variable for A2Plus and A2Minus
HPlus = -kDiff*y + kprolPlus*WNT/(kT*kP + kapPlus);

HMinus = -y*kT*kP/(kprolMinus - kapMinus);

if shouldPlot == 1
    % Plot phase plane results
    figure(); hold on;
    plot(HPlus,y,'k-','linewidth',2);
    plot(y,HMinus,'k-','linewidth',2)
    xlim([0 .2]); ylim([0 2])
    set(gca,'fontsize',18,'linewidth',2); box off
    ylabel('Axin2- Hepatocytes'); xlabel('Axin2+ Hepatocytes')
    title('WNT Dependent Proliferation','fontSize',24)
    
    % Add arrows to plot
%     for A2Plus = 0:.01:.2
%         for A2Minus = 0:.1:2
%                 r1 = A2Plus*kprolPlus*WNT/(A2Plus + kDiff*A2Minus); % A2Plus Hepatocyte proliferation
%                 r2 = A2Plus*kT*kP; % Transition from A2Plus to A2Minus
%                 r3 = A2Plus*kapPlus; % A2Plus Hepatocyte apoptosis
%                 r4 = A2Minus*kprolMinus; % A2Minus Hepatocyte proliferation
%                 r5 = A2Minus*kapMinus; % A2Minus Hepatocyte apoptosis
%                 dxdt(1) = r1 - r2 - r3; % A2+ Hepatocytes
%                 dxdt(2) = r4 + r2 - r5; % A2- Hepatocytes
%                 arrow([A2Plus,A2Minus],[A2Plus+dxdt(1),A2Minus+dxdt(2)])
%         end
%     end
    
    % Plot parameter values vs. no feedback
    figure(); hold on; plot(k3,k0,'ko','markerFaceColor','k','markerSize',4)
    plot([0 1.5],[0 1.5],'k-','linewidth',1); xlim([0 1])
    set(gca,'fontsize',18,'linewidth',2); box off
    xlabel('Implicit competition parameter values')
    ylabel('No feedback parameter values')
    title('WNT Dependent Proliferation','fontsize',24)

end

%% Transition of Axin2+ cells is dependent on WNT concentration (Condition 4)
% Set aliases
A2PlusIC = x0(1); A2MinusIC = x0(2);
kprolPlus = k4(1); kCapPlus = k4(2); kapPlus = k4(6);
kCapMinus = k4(4); kprolMinus = k4(5); kapMinus = k4(7); 
kP = k4(3); kT = k4(8);
WNT = k4(9); kDiff = k4(10);
kRenew = k4(13); C1 = k4(11); C2 = k4(12);

% Equations
y = 0:.01:2; % Dummy variable for A2Plus and A2Minus
HPlus = -kDiff*y + (WNT/(kT*kP))*(kprolPlus - kapPlus);

HMinus = -y.^2*(kT*kP/WNT)./(kprolMinus - kapMinus + y*(kT*kP*kDiff)/WNT);

if shouldPlot == 1
    % Plot phase plane results
    figure(); hold on;
    plot(HPlus,y,'k-','linewidth',2);
    plot(y,HMinus,'k-','linewidth',2)
    xlim([0 .2]); ylim([0 2])
    set(gca,'fontsize',18,'linewidth',2); box off
    ylabel('Axin2- Hepatocytes'); xlabel('Axin2+ Hepatocytes')
    title('WNT Dependent Transition','fontSize',24)
    
    % Add arrows to plot
%     for A2Plus = 0:.01:.2
%         for A2Minus = 0:.1:2
%                 r1 = A2Plus*kprolPlus; % A2Plus Hepatocyte proliferation
%                 r2 = A2Plus*kT*kP*(A2Plus + kDiff*A2Minus)/WNT; % Transition from A2Plus to A2Minus
%                 r3 = A2Plus*kapPlus; % A2Plus Hepatocyte apoptosis
%                 r4 = A2Minus*kprolMinus; % A2Minus Hepatocyte proliferation
%                 r5 = A2Minus*kapMinus; % A2Minus Hepatocyte apoptosis
%                 dxdt(1) = r1 - r2 - r3; % A2+ Hepatocytes
%                 dxdt(2) = r4 + r2 - r5; % A2- Hepatocytes
%                 arrow([A2Plus,A2Minus],[A2Plus+dxdt(1),A2Minus+dxdt(2)])
%         end
%     end
    
    % Plot parameter values vs. no feedback
    figure(); hold on; plot(k3,k0,'ko','markerFaceColor','k','markerSize',4)
    plot([0 1.5],[0 1.5],'k-','linewidth',1); xlim([0 1])
    set(gca,'fontsize',18,'linewidth',2); box off
    xlabel('Implicit competition parameter values')
    ylabel('No feedback parameter values')
    title('WNT Dependent Transition','fontsize',24)

end

%% Alternate population renewing Axin2+ cells (ex. stem cells) (Condition 5)
% Set aliases
A2PlusIC = x0(1); A2MinusIC = x0(2);
kprolPlus = k5(1); kCapPlus = k5(2); kapPlus = k5(6);
kCapMinus = k5(4); kprolMinus = k5(5); kapMinus = k5(7); 
kP = k5(3); kT = k5(8);
WNT = k5(9); kDiff = k5(10);
kRenew = k5(13); C1 = k5(11); C2 = k5(12);

% Equations
y = 0:.1:2; % Dummy variable for A2Plus and A2Minus
HMinus = -y*kT*kP/(kprolMinus - kapMinus);
% Call optimization function to find HPlus
solverfun = @(x)abs(kRenew/(1+exp(C1*x+C2)) + x*(kprolPlus - kT*kP - kapPlus));
[HPlus,fval] = fminsearch(solverfun,1);

if shouldPlot == 1
    % Plot phase plane results
    figure(); hold on;
    plot(HPlus*ones(length(y)),y,'k-','linewidth',2);
    plot(y,HMinus,'k-','linewidth',2)
    xlim([0 .1]); ylim([0 2])
    set(gca,'fontsize',18,'linewidth',2); box off
    ylabel('Axin2- Hepatocytes'); xlabel('Axin2+ Hepatocytes')
    title('Alternate Axin2+ Renewal','fontSize',24)
    
    % Add arrows to plot
%     for A2Plus = 0:.01:.1
%         for A2Minus = 0:.1:2
%                 r1 = A2Plus*kprolPlus + kRenew*(1/(1+exp(C1*A2Plus+C2))); % A2Plus Hepatocyte proliferation
%                 r2 = A2Plus*kT*kP; % Transition from A2Plus to A2Minus
%                 r3 = A2Plus*kapPlus; % A2Plus Hepatocyte apoptosis
%                 r4 = A2Minus*kprolMinus; % A2Minus Hepatocyte proliferation
%                 r5 = A2Minus*kapMinus; % A2Minus Hepatocyte apoptosis
%                 dxdt(1) = r1 - r2 - r3; % A2+ Hepatocytes
%                 dxdt(2) = r4 + r2 - r5; % A2- Hepatocytes
%                 arrow([A2Plus,A2Minus],[A2Plus+dxdt(1),A2Minus+dxdt(2)])
%         end
%     end
    
    % Plot parameter values vs. no feedback
    figure(); hold on; plot(k5,k0,'ko','markerFaceColor','k','markerSize',4)
    plot([0 1.5],[0 1.5],'k-','linewidth',1);
    set(gca,'fontsize',18,'linewidth',2); box off
    xlim([0 1]); ylim([0 1]);
    xlabel('Implicit competition parameter values')
    ylabel('No feedback parameter values')
    title('Alternate Axin2+ Renewal','fontsize',24)

end

end