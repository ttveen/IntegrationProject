%% StateSpace system identification, using parameter fitting
% Add some paths
addpath('..\')
addpath('..\')

% Include tclab.m for initialization
%tclab;

% Load date
clear all
close all
load('..\..\Data\blockResponse5.mat');
time.block = time.long
t1.block = t1.long
t2.block = t2.long
%%

%Ensure that the time vector is evenly spaced
time.sim = linspace(0,time.block(end),length(time.block));
t1.inter = interp1(time.block,t1.block,time.sim);
t2.inter = interp1(time.block,t2.block,time.sim);
u1 = interp1(time.block,u1,time.sim);
u2 = interp1(time.block,u2,time.sim);
h = time.sim(2) - time.sim(1); %Time step

%Input/output corresponding to the data
u = [u1;u2];
y = [t1.inter',t2.inter'];
%Initialisaton
T0 = 273.15; %Celsius to kelvin
Tinf = T0 + 23; %Ambient temperature

T.H1 = zeros(length(time.sim) - 2,1);
T.H2 = zeros(length(time.sim) - 2,1);
T.H1(1) = y(1,1);
T.H2(1) = y(2,1);
T.C1 = y(:,1);
T.C2 = y(:,2);
lsq.y = zeros(2*(length(time.sim) - 2),1);
lsq.F = zeros(2*(length(time.sim) - 2),5);

sigma = 5.67e-8;
tau = 15;
alpha1 = 0.01;
alpha2 = 0.0075;


for i = 1 : length(time.sim) - 2
    %Heater temperauter, expressed in the measured temperature
    T.H1(i+1) = tau/h*(T.C1(i+2)) + (1-tau/h)*(T.C1(i+1));
    T.H2(i+1) = tau/h*(T.C2(i+2)) + (1-tau/h)*(T.C2(i+1));
    
    %Left side of the least square function
    lsq.y(2*i-1,1) = T.H1(i+1) - T.H1(i);
    lsq.y(2*i,1) = T.H2(i+1) - T.H2(i);
    
    %Matrix on the right side of the equation
    lsq.F(2*i-1,:) = h*[Tinf - (T.H1(i)+T0), sigma*(Tinf^4 - (T.H1(i)+T0)^4), (T.H2(i)+T0) - (T.H1(i)+T0), sigma*((T.H2(i)+T0)^4 - (T.H1(i)+T0)^4), alpha1*u1(i)];
    lsq.F(2*i  ,:) = h*[Tinf - (T.H2(i)+T0), sigma*(Tinf^4 - (T.H2(i)+T0)^4), (T.H1(i)+T0) - (T.H2(i)+T0), sigma*((T.H1(i)+T0)^4 - (T.H2(i)+T0)^4), alpha2*u2(i)];
end
lsq.x = (lsq.F'*lsq.F)\(lsq.F'*lsq.y);

%% Compute the model
Tinf = 273.15+23;
TC1 = [20+273.15];
TC2 = [20+273.15];
f.y = [20;20;20;20]+ 273.15;
t = 0;
for i = 1: length(time.sim)
   f = ode45(@(t,x)heaters(t,x,u1(i),u2(i),lsq, alpha1, alpha2, tau),[0 h],f.y(:,end));
   TC1 = [TC1, f.y(3,2:end)];
   TC2 = [TC2, f.y(4,2:end)];
   %t = [t, f.x(2:end)];
end
TC1 = TC1 - 273.15;
TC2 = TC2 - 273.15;
t = linspace(0,time.block(end), length(TC1));
figure(stepBlock)
subplot(2,1,1)
hold on
plot(t,TC1)
plot(t,TC2)
legend({'Heater 1', 'Heater 2','H1, lsq fit','H1, lsq fit'},'Location','northeast')
saveas(stepBlock,'..\..\latex\images\SYSID\leastSquaresFit','svg')


%% save workspace variables
save('../../Data/LeastSquares.mat', 't', 'TC1', 'TC2', 'lsq');

%% function 
function dx = heaters(t,x,u1,u2,lsq, alpha1, alpha2, tau)
    sigma = 5.67e-8;
    %tau = 23;
    %alpha1 = 0.01;
    %alpha2 = 0.0075;
    T0 = 273.15; %Celsius to kelvin
    Tinf = T0 + 23; %Ambient temperature

    dx(1,1) = lsq.x(1)*(Tinf - x(1)) + lsq.x(2)*sigma*(Tinf^4 - x(1)^4) + lsq.x(3)*(x(2)-x(1)) + sigma*lsq.x(4)*(x(2)^4 - x(1)^4) + alpha1*u1;
    dx(2,1) = lsq.x(1)*(Tinf - x(2)) + lsq.x(2)*sigma*(Tinf^4 - x(2)^4) - lsq.x(3)*(x(2)-x(1)) - sigma*lsq.x(4)*(x(2)^4 - x(1)^4) + alpha2*u2;
    dx(3,1) = (x(1) - x(3))/tau;
    dx(4,1) = (x(2) - x(4))/tau;
end
