%% MPC Controller Design
clear; clc; close all;


%% Load and initialise experiment data
addpath('../')
load('../../Data/NSID');
tclab;

A = NSID.At;
B = NSID.Bt;
C = NSID.Ct;
D = NSID.Dt;
K = NSID.Ks;
%% Init variables
% init Time
runtime = 1000; %1000 second runtime
% axis([fig1a fig1b], [0 runtime 0 50]);
t = zeros(1,runtime);
MPCcompTime = t;
time.MPC = zeros(1,runtime);

% init Kalman filter
xhat = zeros(size(A,1),runtime);
yhat = zeros(2,runtime);

% init measurements
t1.MPC = zeros(1,runtime);
t2.MPC = zeros(1,runtime);

% init input
u = zeros(2,runtime);
%% MPC

% MPC parameters
N = 80; % Horizon
nu = 2;
nx = size(A,1);
ny = size(C,1);

% input and output bounds
ulb = 0;
uub = 50;
ylb = 0;
yub = 60;

%%
Q = eye(2);
R = 0.001*eye(2);
P = eye(2);
% Create symbolic decision variables
u_ = sdpvar(repmat(nu,1,N),repmat(1,1,N));
x = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));
y = sdpvar(repmat(ny,1,N),repmat(1,1,N));
r = sdpvar(repmat(ny,1,N+1),repmat(1,1,N+1));

% initialisation of constraints and objective variables
constraints = [];
objective = 0;

for k = 1:N-1
    % objective function, Quadratic over the state and input
    objective = objective + (y{k}-r{k})'*Q*(y{k}-r{k}) + (u_{k})'*R*(u_{k});
    
    % constraints, first the dynamics, then the input output bounds
    constraints = [constraints, x{k+1} == A*x{k}+B*u_{k}];
    constraints = [constraints, y{k} == C*x{k}+D*u_{k}];
    constraints = [constraints, ulb <= u_{k} <= uub];
    constraints = [constraints, ylb <= y{k} <= yub];
end

% The terminal cost
objective = objective + (y{N}-r{N})'*P*(y{N}-r{N});

parameters_in = {x{1},[r{:}]};
solutions_out = u_{:};
%The controller, given the constraints and objective
controller = optimizer(constraints,objective,sdpsettings('solver','gurobi'),parameters_in,solutions_out);

x = zeros(6,1);
MPC_x(:,1) = x;
MPC_u(:,1) = [0;0];

for i = 1:runtime-1
    future_r = [35+3*sin((i:i+N)/100); 30+3*cos((i:i+N)/150)];
    inputs = {x,future_r};
    U = controller{inputs};
%     U = solutions{1};
    %X = solutions{2};
    x = A*x+B*U(:,1);
   
   MPC_u(:,i+1) = U(:,1);
   MPC_x(:,i+1) = x;
end

y = C*MPC_x+D*MPC_u;
yr = [35+3*sin((0:1:runtime)/100); 30+3*cos((0:1:runtime)/150)];
figure
hold on
plot(yr(1,:),'r--')
plot(yr(2,:),'b--')
plot(y(1,:))
plot(y(2,:))
hold off

figure
hold on
plot(MPC_u(1,:))
plot(MPC_u(2,:))
hold off
%% Experiment
figure('Name','MPCexp')
fig1a = subplot(2,1,1);
aniT1 = animatedline(fig1a,'Color','r','Marker','.','Linestyle','none');
aniT2 = animatedline(fig1a,'Color','b','Marker','.','Linestyle','none');
aniK1 = animatedline(fig1a,'Color','r','Marker','+','Linestyle','none');
aniK2 = animatedline(fig1a,'Color','b','Marker','+','Linestyle','none');
hold on
xlabel('Time in s')
ylabel({'Temperature in $^{\circ}C$'},  'Interpreter', 'Latex')
title('Temperature')
axis([0 runtime 0 45])
fig1b = subplot(2,1,2);
aniu1 = animatedline(fig1b,'Color','r','Marker','.','Linestyle','none');
aniu2 = animatedline(fig1b,'Color','b','Marker','.','Linestyle','none');
xlabel('Time in s')
ylabel('Input level',  'Interpreter', 'Latex')
title('Input')
axis([0 runtime -5 50])
hold on
%Reference signal
yr = [35+3*sin((0:1:runtime)/100); 30+3*cos((0:1:runtime)/150)];
plot(fig1a,(0:1:runtime),yr(1,:),'r--')
plot(fig1a,(0:1:runtime),yr(2,:),'b--')
led(0.5)
for i = 1:runtime
    currenttime = clock;
    tic
    % generate reference
    future_r = [35+3*sin((i:i+N)/100); 30+3*cos((i:i+N)/150)];
    % generated input
    tic
    inputs = {xhat(:,i),future_r};
    U = controller{inputs};
    u(:,i) = U(:,1);
    % implement input
    h1(u(1,i));
    h2(u(2,i));
    
    % measure output
    t1.MPC(i) = T1C();
    t2.MPC(i) = T2C();
    
    % Kalman filter
    xhat(:,i+1) = (A-K*C)*xhat(:,i) + B*[u(1,i); u(2,i)] + K*[t1.MPC(i); t2.MPC(i)];
    yhat(:,i) = C*xhat(:,i) + D*[u(1,i); u(2,i)];
    MPCcompTime(i) = toc; %Measure what for time it takes to loop through the whole system
    
    % plot results
    addpoints(aniT1,time.MPC(i),t1.MPC(i))
    %plot(fig1a,time.MPC(i),t1.MPC(i),'r.')
    hold on
    addpoints(aniT2,time.MPC(i),t2.MPC(i))
    addpoints(aniK1,time.MPC(i),yhat(1,i))
    addpoints(aniK2,time.MPC(i),yhat(2,i))
    
    addpoints(aniu1, time.MPC(i), u(1,i))
    hold on
    addpoints(aniu2, time.MPC(i), u(2,i))
    
    % time stuff
    t(i) = toc;
    pause(max(0.01,1-t(i)));
    time.MPC(i+1) = time.MPC(i) + etime(clock,currenttime);
end

turnOff(a)
led(0)
time.MPC(end) =[];
yr(:,end) = [];


figure('Name','MPC, peridiodic reference')
fig2a = subplot(2,1,1);
plot(fig2a, time.MPC, t1.MPC, 'r.','MarkerSize',8)
hold on
plot(fig2a, time.MPC, t2.MPC, 'b.','MarkerSize',8)
plot(fig2a, time.MPC, yhat(1,:), 'r','LineWidth',3)
plot(fig2a, time.MPC, yhat(2,:), 'b','LineWidth',3)
plot(fig2a, time.MPC, yr(1,:), 'r--','LineWidth',3)
plot(fig2a, time.MPC, yr(2,:), 'b--','LineWidth',3)
xlabel('Time in s')
ylabel({'Temperature in $^{\circ}C$'},  'Interpreter', 'Latex')
title('Temperature')
axis(fig2a,[0 runtime -inf inf])
legend(fig2a,{'$T_{C1,measured}$','$T_{C2,measured}$','$T_{C1,estimate}$','$T_{C2,estimate}$','$T_{C1,reference}$','$T_{C1,reference}$'},'Interpreter','latex','Location','southeast','FontSize', 20)
fig2b = subplot(2,1,2);
xlabel('Time in s')
ylabel('Input level',  'Interpreter', 'Latex')
title('Input')
hold on
plot(fig2b, time.MPC, u(1,:), 'r.','MarkerSize',8)
plot(fig2b, time.MPC, u(2,:), 'b.','MarkerSize',8)
legend(fig2b,{'$u_1$', '$u_2$'}, 'Interpreter','latex','FontSize', 20)
axis(fig2b,[0 runtime -inf inf])

MPCexpper1 = gcf;
set(MPCexpper1, 'position', get(0, 'ScreenSize'))
MPCexpper1.Renderer = 'painters';
saveas(MPCexpper1, '../../Latex/images/controller/MPCper1', 'svg');
save('../../Data/MPCper1.mat', 'time', 't1', 't2', 'MPCexpper1', 'u','xhat','yhat', 'N','MPCcompTime')
close gcf