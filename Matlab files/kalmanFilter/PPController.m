%% Kalman filter test

%% Load and initialise experiment data
addpath('../')
load('../../Data/NSID');
tclab;
A = NSID.At;
B = NSID.Bt;
C = NSID.Ct;
D = NSID.Dt;
K = NSID.Ks;

%% Create figure
figure('Name','Kalman Filter Test')
fig1a = subplot(2,1,1);
hold on
xlabel('Time in s')
ylabel({'Temperature in $^{\circ}C$'},  'Interpreter', 'Latex')
title('Temperature')
fig1b = subplot(2,1,2);
xlabel('Time in s')
ylabel('Input level',  'Interpreter', 'Latex')
title('Input')
hold on

%% Init variables
% init Time
runtime = 400; %200 second runtime
% axis([fig1a fig1b], [0 runtime 0 50]);
t = zeros(1,runtime);
time.kalman = zeros(1,runtime);

% init Kalman filter
xhat = zeros(size(A,1),runtime);
yhat = zeros(2,runtime);

% init measurements
t1.kalman = zeros(1,runtime);
t2.kalman = zeros(1,runtime);

% init input
u = zeros(2,runtime);

% init reference
r = 50*ones(2,runtime);

%% Determine gain
% p = [-0.9 -0.8 -0.7 0.7 0.8 0.9]; % stable discrete poles
% p = [-0.3 -0.2 -0.1 0.1 0.2 0.3];
p = 0.2*[0.7 0.75 0.8 0.85 0.9 0.95];

Kx = place(A,B,p); % pole placement

Kr = inv(C/(eye(size(A))-A+B*Kx)*B); % reference gain

%% simulate
x = xhat;
x(:,1) = C\[23; 23]; 
for i = 1:runtime-1
    % generated input
    u(:,i) = -Kx*x(:,i) + Kr*r(:,i);
    
    % system
    x(:,i+1) = A*x(:,i)+B*u(:,i);
        
end
y = C*x+D*u;

%% Run Experiment

%Measure the corresponding output
led(1);

for i = 1:runtime - 1
    currenttime = clock;
    tic
    
    % generated input
    u(:,i) = -Kx*xhat(:,k) + Kr\r(:,k);
    
    % implement input
    h1(u(1,i));
    h2(u(2,i));
    
    % measure output
    t1.kalman(i) = T1C();
    t2.kalman(i) = T2C();
    
    % Kalman filter
    xhat(:,i+1) = (A-K*C)*xhat(:,i) + B*[u1(i); u2(i)] + K*[t1.kalman(i); t2.kalman(i)];
    yhat(:,i) = C*xhat(:,i) + D*[u1(i); u2(i)];
    
    
    % plot results
    plot(fig1a,time.kalman(i),t1.kalman(i),'r.')
    hold on
    plot(fig1a,time.kalman(i),t2.kalman(i),'b.')
    plot(fig1a,time.kalman(i), yhat(1,i), 'r+')
    plot(fig1a,time.kalman(i), yhat(2,i), 'b+')
    %legend(fig2a,'Heater 1','Location','northwest')
    
    plot(fig1b, time.kalman(i), u1(i), 'r.')
    hold on
    plot(fig1b, time.kalman(i), u2(i), 'b.')
    %legend(fig2b,'Heater 2','Location','northwest')
    
    % time stuff
    t(i) = toc;
    pause(max(0.01,1-t(i)));
    time.kalman(i+1) = time.kalman(i) + etime(clock,currenttime);
end

turnOff(a)
led(0)
time.kalman(end) = [];

kalmanTestFig = gcf;
stepBlock.Renderer = 'painters';
saveas(kalmanTestFig, '../../Latex/images/kalmanTest/kalmanTest2', 'svg');
save('../../Data/kalmanTest2.mat', 'time', 't1', 't2', 'kalmanTestFig', 'u1', 'u2', 'xhat', 'yhat')