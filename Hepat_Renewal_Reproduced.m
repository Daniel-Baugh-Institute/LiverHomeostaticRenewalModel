function [t,cell] = FileS4_HepatRenewalReproduced()
% global kdhigh_vec kdlow_vec trange h_iter time_iter iter
global kdhigh kdlow kdhigh_vec kdlow_vec trange h iter time
iter=1;

k=zeros(17,1);
%Table A.1. Parameter values
k(1) = 0.2143;  %kprolhigh
k(2) = 0.075;   %kcaphigh   
k(3) = 0.435;   %kenvprol
k(4) = 0.1;     %kA
k(5) = 0.0446;  %kdeathhigh
k(6) = 0.17;    %kT
k(7) = 0.95; %0.95 in Dan's code and is recalculated; %0.17 in Table A1;    %kenvT
k(8) = 8;       %krenew
k(9) = 100;     %C1
k(10) = 8;      %C2
k(11) = 0.1071; %kprollow
k(12) = 1.425;  %kcaplow
k(13) = 0.0446; %kdeathlow
k(14) = 37.5;   %kp1
k(15) = 20.8;   %ki1
k(16) = 15.9;   %kp2
k(17) = 2.5;    %ki2
%Table ends here

kdhigh=k(5);
kdlow=k(13);

%define time range
tStart = 0;
tEnd = 100;

trange=tStart:0.1:tEnd;
kdhigh_vec=k(5)*ones(1,length(trange));
kdlow_vec=k(13)*ones(1,length(trange));

kdhigh_vec(floor(mod(trange,2))==1)=2*k(5);
kdlow_vec(floor(mod(trange,2))==1)=2*k(13);

%set the initial steady state condition for solving the model equations
cell0 = [0.05,0.95,0,0];

%k
%RenewalODE(0,[0.05,0.95,0,0],k,[0.05,0.95,0,0])

%calling ode solver to solve the system of equations
[t,cell] = ode15s(@(t,cell)RenewalODE(t,cell,k,cell0), [tStart tEnd], cell0);

save('var.mat');

% figure;box on;
% plot()

figure;box on;
plot(t,cell(:,1));
xlabel('Time (days)');
ylabel('SR_{high} Cells');

figure;box on;
plot(t,cell(:,2))
xlabel('Time (days)');
ylabel('SR_{low} Cells');

figure;box on;
plot(cell(:,1),cell(:,2))
xlabel('SR_{high} Cells');
ylabel('SR_{low} Cells');

figure;box on;
plot(t,[cell(:,3),cell(:,4)])
xlabel('Time (days)');
ylabel('Control Action');


end

%%
function [dxdt] = RenewalODE(t,x,k,cell0)

% global kdhigh_vec kdlow_vec trange iter time_iter h_iter
% dcelldt=zeros(2,1);

global h iter kdhigh kdlow trange kdhigh_vec kdlow_vec time

SRhigh = x(1); SRlow = x(2); 
CA1 = x(3); CA2 = x(4);
% CA1=0;CA2=0;
SRhigh0 = 0.05;SRlow0 = 0.95;

% kdhigh=k(5);
% kdlow=k(13);
% kdhigh=kdhigh_vec(find(trange==floor(t),1))
% kdlow=kdlow_vec(find(trange==floor(t),1))

% time_iter(iter)=t;


% kdhigh(floor(mod(t,2))==1)=2*k(5);
% kdlow(floor(mod(t,2))==1)=2*k(13);
% 
% k(5)=kdhigh;k(13)=kdlow;

tp=find(t==trange,1);


% if ~isempty(tp)
%     k(5)=kdhigh_vec(tp);
%     k(13)=kdhigh_vec(tp);
% 
% end

% if floor(mod(t,2))==1
%     k(5)=2*k(5);
%     k(13)=2*k(13);
% end

rem=mod(t,2);
if rem >= 1
    k(5)=1.5*kdhigh;
    k(13)=1.5*kdlow;
else
    k(5)=kdhigh;
    k(13)=kdlow;
    kpseudo=k(5);
end

% h(iter)=k(5);
% iter=iter+1;

if ~isempty(tp)
    time(iter)=t;
    h(iter)=k(5);
    iter=iter+1;
end


% k

% Calculate reaction rates
r1 = SRhigh * k(1) * (1-(SRhigh/k(2))) * ((k(3)*(1+CA1))/(SRhigh+k(4)*SRlow));
r2 = SRhigh * k(5);
r3 = SRhigh * k(6) * ((k(7)*(1+CA2))/(SRlow));
r4 = k(8) * (1/(1+exp(k(9)*SRhigh+k(10))));
r5 = SRlow * k(11) * (1-(SRlow/k(12)));
r6 = SRlow * k(13); 

r7 = k(15) * (SRhigh0 - SRhigh);
r8 = k(17) * (SRlow0 - SRlow);

%defining the  model equations
dxdt(1) = r1 - r2 - r3 + r4;
dxdt(2) = r3 + r5 - r6;

dxdt(3) = -k(14) * dxdt(1) + r7;
dxdt(4) = -k(16) * dxdt(2) + r8;

dxdt=dxdt';
end