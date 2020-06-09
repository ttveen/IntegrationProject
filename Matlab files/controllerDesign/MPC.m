%% MPC Controller Design
clear; clc; close all;


%% Load and initialise experiment data
addpath('../')
load('../../Data/NSID');
%tclab;

A = NSID.At;
B = NSID.Bt;
C = NSID.Ct;
D = NSID.Dt;
K = NSID.Ks;

%% Init variables
% init Time
runtime = 400; % seconds

%% MPC

% MPC parameters
N = 100; % Horizon
nu = 2;
nx = size(A,1);
ny = size(C,1);

% input and output bounds
ulb = 0;
uub = 50;
ylb = 0;
yub = 60;

% solution to the discrete time Riccati equations
% Q = C'*C;
R = diag([0.3 0.3]);
% % R = 0.3;
% 
% [P,F,sys_p] = idare(A,B,Q,R);
% 
% I = eye(size(A));
% G = inv((C-D*F)*inv(I-A+B*F)*B+D);

%% Reference stuff
N = runtime;
r = 30*ones(2,N);
cvx_begin quiet
    variables x_r(6,N+1) u_r(2,N) y_r(2,N)
    minimize sum(norm(y_r-r,2))
    subject to
%     y_r(:,1) == 23;
    for k = 1:N
       x_r(:,k+1) == A*x_r(:,k) + B*u_r(:,k);
       y_r(:,k) == C*x_r(:,k) + D*u_r(:,k);
       0 <= u_r(:,k) <= 50;      
    end
cvx_end

%%
Q = eye(2);
R = eye(2);
P = eye(2);
% Create symbolic decision variables
u = sdpvar(repmat(nu,1,N),repmat(1,1,N));
x = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));
y = sdpvar(repmat(ny,1,N),repmat(1,1,N));

% initialisation of constraints and objective variables
constraints = [];
objective = 0;
% 
% Q1 = 1000*eye(2);
% R1 = eye(2);

% P = idare(A,B,C'*Q1*C,D'*R1*D);
% P = 1000*eye(2);
% r = [20; 20];
for k = 1:N-1
    % objective function, Quadratic over the state and input
    objective = objective + (y{k}-[30;30])'*Q*(y{k}-[30;30]) + (u{k})'*R*(u{k});
%     objective = objective + (y{k}-r)'*Q1*(y{k}-r) + u{k}'*R*u{k};
%     objective = objective + x{k}'*C'*Q1*C*x{k} + u{k}'*(D'*Q1*D+R1)*u{k}...
%         + x{k}'*C'*Q1*D*u{k} + u{k}'*D'*Q1*C*x{k} ...
%         -x{k}'*C'*Q1*r(:,1)... 
%         + u{k}'*D'*Q1*r(:,1)...
%         -r(:,1)'*Q1*C*x{k} ...
%         + r(:,1)'*Q1*r(:,1) + u{k}'*R1*u{k};
%     objective = objective + norm(y{k}-r(:,1),2) + norm(u{k},2);
    
    % constraints, first the dynamics, then the input output bounds
%     constraints = [constraints, xhat{k+1}==(A-K*C)*xhat{k}+B*u{k}+K*y{k}];
    constraints = [constraints, x{k+1} == A*x{k}+B*u{k}];
    constraints = [constraints, y{k} == C*x{k}+D*u{k}];
    constraints = [constraints, ulb <= u{k} <= uub];
    constraints = [constraints, ylb <= y{k} <= yub];
end

% The terminal cost
objective = objective + (y{N}-[30;30])'*P*(y{N}-[30;30]);
% objective = objective + x{N}'*C'*P*C*x{N}...
%     + u{N}'*D'*P*C*x{N}...
%     + x{N}'*C'*P*D*u{N} + u{N}'*D'*P*D*u{N} - x{N}'*C'*P*r(:,1)...
%     - u{N}'*P'*P*r(:,1) - r(:,1)'*P*C*x{N} - r(:,1)'*P*D*u{N} + r(:,1)'*P*r(:,1);

%The controller, given the constraints and objective
controller = optimizer(constraints,objective,[],x{1},[u{:}]);

y0 = [23; 23];
% x = (C-D*F)\(y0-D*G*r);
x = zeros(6,1);
MPC_x(:,1) = x;
MPC_u(:,1) = [0;0];

for i = 1:runtime-1
    U = controller{x};
    x = A*x+B*U(:,1);
   
   MPC_u(:,i+1) = U(:,1);
   MPC_x(:,i+1) = x;
   i
end

y = C*MPC_x+D*MPC_u;

% figure
% hold on
% plot(MPC_x(1,:))
% plot(MPC_x(2,:))
% plot(MPC_x(3,:))
% plot(MPC_x(4,:))
% plot(MPC_x(5,:))
% plot(MPC_x(6,:))
% hold off

figure
hold on
plot(y(1,:))
plot(y(2,:))
hold off

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


