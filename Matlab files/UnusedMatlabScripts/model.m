clear; clc;

% Parameters
parm.T0 = [296.15, 296.15];  
parm.Tinf = 23 + 273.15;   
parm.U = 10.0;           
parm.m = 4.0/1000.0;     
parm.Cp = 0.5 * 1000.0;  
parm.A = 10.0 / 100.0^2; 
parm.As = 2.0 / 100.0^2; 
parm.a1 = 0.0100;    
parm.a2 = 0.0075;    
parm.eps = 0.9;          
parm.sigma = 5.67e-8;    

% simulation time
T = 600;
dt = 0.5;
time = 0:dt:T; % Time vector
tn = numel(time);

% heater input [0 100]%
Q1 = ones(tn,1)*50;
Q2 = ones(tn,1)*50;

% initialise temperature vectors
T1 = ones(tn,1)*parm.T0(1);
T2 = ones(tn,1)*parm.T0(2);

for i = 2:tn
    x0 = [T1(i-1),T2(i-1)];
    tStep = [time(i-1),time(i)];
    
    % Simulate system for one time step
    f = ode45(@(t,x)heater(t,x,Q1(i-1),Q2(i-1),parm),tStep,x0);
    
    T1(i) = f.y(1,end);
    T2(i) = f.y(2,end);
end

% Convert Temperatures from Kelvin to Celcius
T1c = T1 - 273.15;
T2c = T2 - 273.15;

% heater model
function dTdt = heater(t,x,Q1,Q2,parm)

    QC12 = parm.U*parm.As*(x(2)-x(1));
    QR12 = parm.eps*parm.sigma*parm.A*(x(2)^4-x(1)^4);
    dT1dt = (1/(parm.m*parm.Cp))*parm.U*parm.A*(parm.Tinf-x(1))...
        + parm.eps*parm.sigma*parm.A*(parm.Tinf^4-x(1)^4)...
        + QC12 + QR12 + parm.a1*Q1;
    dT2dt = (1/(parm.m*parm.Cp))*parm.U*parm.A*(parm.Tinf-x(2))...
        + parm.eps*parm.sigma*parm.A*(parm.Tinf^4-x(2)^4)...
        - QC12 - QR12 + parm.a2*Q2;
    
    dTdt = [dT1dt, dT2dt]';
end