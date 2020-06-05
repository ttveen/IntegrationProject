%% StateSpace system identification, using parameter fitting
% Add some paths
addpath('..\')
addpath('..\')

% Include tclab.m for initialization
%tclab;

% Load date
clear all
close all
load('..\..\Data\blockResponse2.mat');
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
tau = 10;
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
%%
load('..\..\Data\blockResponse.mat');
%Ensure that the time vector is evenly spaced
time.h1 = linspace(0,325,326);
time.sim = linspace(0,time.block(end),length(time.h1));
t1.inter2 = interp1(time.block,t1.block,time.h1);
t2.inter2 = interp1(time.block,t2.block,time.h1);
u1 = interp1(time.block,u1,time.h1);
u2 = interp1(time.block,u2,time.h1);
h = time.h1(2) - time.h1(1); %Time step

%Input/output corresponding to the data
u = [u1;u2];
y = [t1.inter2',t2.inter2'];
%
T0 = 273.15;
sigma = 5.67e-8;
tau = 23;
alpha1 = 0.01;
alpha2 = 0.0075;
SimTime = 326;
Ta = 273+23;
Temperature = zeros(SimTime,3);
Temperature(:,2) = y(:,1);
Temperature(:,3) = y(:,2);
y1 = zeros(SimTime-2,1);
y2 = zeros(SimTime-2,1);
y = zeros(2*(SimTime-2),1);
F1 = zeros(SimTime-2,5);
F2 = zeros(SimTime-2,5);
for i= 1:SimTime-2
    y1(i) = tau*Temperature(i+2,2)-(tau-1)*Temperature(i+1,2)-(tau*Temperature(i+1,2)-(tau-1)*Temperature(i,2));
    y2(i) = tau*Temperature(i+2,3)-(tau-1)*Temperature(i+1,3)-(tau*Temperature(i+1,3)-(tau-1)*Temperature(i,3));
    %y(2*i-1:2*i) = [y1(i);y2(i)];
    F1(i,:) = [Ta-(tau*(Temperature(i+1,2)+T0)-(tau-1)*(Temperature(i,2)+T0)) sigma*(Ta^4-(tau*(Temperature(i+1,2)+T0)-(tau-1)*(Temperature(i,2)+T0))^4) ...
        (tau*(Temperature(i+1,3)+T0)-(tau-1)*(Temperature(i,3)+T0))-(tau*(Temperature(i+1,2)+T0)-(tau-1)*(Temperature(i,2)+T0)) sigma*((tau*(Temperature(i+1,3)+T0)-(tau-1)*(Temperature(i,3)+T0))^4 - (tau*(Temperature(i+1,2)+T0)-(tau-1)*(Temperature(i,2)+T0))^4) ...
        alpha1*u1(i)];
    F2(i,:) = [Ta-(tau*(Temperature(i+1,3)+T0)-(tau-1)*(Temperature(i,3)+T0)) sigma*(Ta^4-(tau*(Temperature(i+1,3)+T0)-(tau-1)*(Temperature(i,3)+T0))^4) ...
        (tau*(Temperature(i+1,2)+T0)-(tau-1)*(Temperature(i,2)+T0))-(tau*(Temperature(i+1,3)+T0)-(tau-1)*(Temperature(i,3)+T0)) sigma*((tau*(Temperature(i+1,2)+T0)-(tau-1)*(Temperature(i,2)+T0))^4 - (tau*(Temperature(i+1,3)+T0)-(tau-1)*(Temperature(i,3)+T0))^4) ...
        alpha2*u2(i)];
end
y = [y1;y2];
F = [F1;F2];
lsq.x = inv(F'*F)*F'*y;
%%
lsq.x= [5.53*0.001/0.5; 0.9*0.001/0.5; 24.44*0.0002/0.5; 0.9*0.0002/0.5; 1/0.05];
%% Compute the model
Tinf = 273.15+23;
TC1 = [298];
TC2 = [298];
f.y = [298;298;298;298];
t = 0;
for i = 1: length(time.sim)
   f = ode45(@(t,x)heaters(t,x,u1(i),u2(i),lsq),[0 h],f.y(:,end));
   TC1 = [TC1, f.y(3,2:end)];
   TC2 = [TC2, f.y(4,2:end)];
   %t = [t, f.x(2:end)];
end
TC1 = TC1 - 273.15;
TC2 = TC2 - 273.15;
figure()
plot(TC1)
hold on
plot(TC2)
function dx = heaters(t,x,u1,u2,lsq)
    sigma = 5.67e-8;
    tau = 23.16;
    alpha1 = 0.01;
    alpha2 = 0.0075;
    T0 = 273.15; %Celsius to kelvin
    Tinf = T0 + 23; %Ambient temperature

    dx(1,1) = lsq.x(1)*(Tinf - x(1)) + lsq.x(2)*sigma*(Tinf^4 - x(1)^4) + lsq.x(3)*(x(2)-x(1)) + sigma*lsq.x(4)*(x(2)^4 - x(1)^4) + alpha1*u1;
    dx(2,1) = lsq.x(1)*(Tinf - x(2)) + lsq.x(2)*sigma*(Tinf^4 - x(2)^4) - lsq.x(3)*(x(2)-x(1)) - sigma*lsq.x(4)*(x(2)^4 - x(1)^4) + alpha2*u2;
    dx(3,1) = (x(1) - x(3))/tau;
    dx(4,1) = (x(2) - x(4))/tau;
end
