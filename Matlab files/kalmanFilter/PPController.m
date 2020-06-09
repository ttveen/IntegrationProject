% MPC Controller Design
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
runtime = 6000; % seconds

% init reference
r = 30*ones(2,runtime);
% r = (20*sin(0.01*(0:runtime-1))+30) .* [1;-1]; 
% r1 = 30+5*sin(0.002*(0:1:runtime));
% r2 = 30+5*cos(0.002*(0:1:runtime));
% r = [r1;r2];

%% MPC

% MPC parameters
N = 20; % Horizon
nu = 2;
nx = size(A,1);
ny = size(C,1);

% input and output bounds
ulb = 0;
uub = 100;
ylb = 0;
yub = 60;

% solution to the discrete time Riccati equations
Q = C'*C;
% R = diag([0.3 0.5]);
R = 0.3;

[P,F,sys_p] = idare(A,B,Q,R);

I = eye(size(A));
G = inv((C-D*F)*inv(I-A+B*F)*B+D);

% Create symbolic decision variables
u = sdpvar(repmat(nu,1,N),repmat(1,1,N));
xhat = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));
y = sdpvar(repmat(ny,1,N),repmat(1,1,N));

% initialisation of constraints and objective variables
constraints = [];
objective = 0;

for k = 1:N-1
    % objective function, Quadratic over the state and input
    objective = objective + xhat{k}'*Q*xhat{k} + u{k}'*R*u{k};
    
    % constraints, first the dynamics, then the input output bounds
    constraints = [constraints, xhat{k+1}==(A-K*C)*xhat{k}+B*u{k}+K*y{k}];
    constraints = [constraints, ulb <= u{k} <= uub];
    constraints = [constraints, ylb <= y{k} <= yub];
end

% The terminal cost
objective = objective + (xhat{N}'*P*xhat{N});

%The controller, given the constraints and objective
controller = optimizer(constraints,objective,[],xhat{1},[u{:}]);

y0 = [23; 23];
xhat = (C-D*F)\(y0-D*G*r(:,1));
MPC_xhat(:,1) = xhat;
MPC_y(:,1) = 0;
MPC_u(:,1) = 0;

for i = 1:runtime-1
   U = controller{xhat};
   y
   xhat = (A-K*C)*xhat+B*U(1)+K*y{k};
end
%% simulate
% x = xhat;
% y0 = [23; 23];
% x(:,1) = (C-D*F)\(y0-D*Kr*r(:,1));
% for i = 1:runtime-1
%     % generated input
%     u(:,i) = -F*x(:,i)+Kr*r(:,i);
%     
%     % system
%     x(:,i+1) = A*x(:,i)+B*u(:,i);
%         
% end
% y = C*x + D*u;
% Y = C*inv(I-A)*B*u(:,end-1)+D*u(:,end-1); %steady-state values
% % figure;
% % plot(y(1,:))
% % hold on
% % plot(y(2,:))
% % hold off
% % 
% % figure;
% % plot(u(1,:))
% % hold on
% % plot(u(2,:))
% % hold off
% 
% plot(y(1,:))
% hold on
% plot(y(2,:))
% plot(r(1,:))
% plot(r(2,:))
% legend('y1','y2','r1','r2')
% figure()
% plot(u(1,:))
% hold on
% plot(u(2,:))
% legend('u1','u2')


