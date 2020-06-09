%% MPC Controller Design

% clear; clc; close all;


%% Load and initialise experiment data

addpath('../')
load('../../Data/NSID');
tclab;

A = NSID.At;
B = NSID.Bt;
C = NSID.Ct;
D = NSID.Dt;
K = NSID.Ks;


%% MPC preparations 

% parameters
nu = 2;
nx = size(A,1);
ny = size(C,1);
N = 80; % Horizon
runtime = 1000; % seconds
% r = [30;30];

% solution to the discrete time Riccati equations
Q = C'*C;
R = 0.01*diag([0.3 0.3]);

[P,F,sys_p] = idare(A,B,Q,R);


%% Determine reachable reference trajectory (step input)

% Nr = N;
% t = 0:1:Nr-1;
% r = ones(2,1).*(1.25*(sin((pi/(Nr/2))*t-0.5*pi)+1))+23;
% cvx_begin quiet
%     variables x_r(6,Nr+1) u_r(2,Nr) y_r(2,Nr)
%     minimize sum(norm(y_r-r,2))
%     subject to
%     for k = 1:Nr
%        x_r(:,k+1) == A*x_r(:,k) + B*u_r(:,k);
%        y_r(:,k) == C*x_r(:,k) + D*u_r(:,k);
%        0 <= u_r(:,k) <= 50;
% %        0 <= y_r(:,k);
%     end
% cvx_end


%% MPC Objective function and constraints initialisation

% % input and output bounds
% ulb = 0;
% uub = 50;
% ylb = 0;
% yub = 60;
% 
% % Create symbolic decision variables
% u = sdpvar(repmat(nu,1,N),repmat(1,1,N));
% x = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));
% y = sdpvar(repmat(ny,1,N),repmat(1,1,N));
% 
% % initialisation of constraints and objective variables
% constraints = [];
% objective = 0;
% 
% for k = 1:N-1
%     % objective function, Quadratic over the state and input
%     objective = objective + (x{k}-x_r(:,k))'*Q*(x{k}-x_r(:,k)) + (u{k}-u_r(:,k))'*R*(u{k}-u_r(:,k));
%     
%     % constraints, first the dynamics, then the input output bounds
%     constraints = [constraints, x{k+1} == A*x{k}+B*u{k}];
%     constraints = [constraints, y{k} == C*x{k}+D*u{k}];
%     constraints = [constraints, ulb <= u{k} <= uub];
%     constraints = [constraints, ylb <= y{k} <= yub];
% end
% 
% % The terminal cost
% objective = objective + (x{N}-x_r(:,N))'*P*(x{N}-x_r(:,N));

%% Poging 2
Nr = 100;
t = 0:1:Nr-1;
r = ones(2,1).*(1.25*(sin((pi/(Nr/2))*t-0.5*pi)+1))+23;
Q = eye(2);
R = 0.01*eye(2);
cvx_begin quiet
    variables x(6,N+1) u(2,N) y(2,N) x_s(6,Nr+1) u_s(2,Nr) y_s(2,Nr)
%     minimize sum((y-y_s(:,1:N))'*Q*(y-y_s(:,1:N)) + (u-u_s(:,1:N))'*R*(u-u_s(:,1:N)))+ sum(norm(y_s-r,2));
    minimize sum(norm(y-y_s(:,1:N),2) + norm(u-u_s(:,1:N),2)) + sum(norm(y_s-r,2));
    subject to
    for k = 1:N
        y(:,k) == C*x(:,k)+D*u(:,k);
        0 <= u(:,k) <= 50;
        0 <= y(:,k) <= 60;
        
        x_s(:,k+1) == A*x_s(:,k) + B*u_s(:,k);
        y_s(:,k) == C*x_s(:,k) + D*u_s(:,k);
        0 <= u_s(:,k) <= 50;
        0 <= y_s(:,k) <= 60;
    end
%     x_s(:,1) == A*x_s(:,Nr-1)+B*u_s(Nr-1);
    x(:,N) == x_s(:,N);
    
cvx_end

%% MPC simulation

% %The controller, given the constraints and objective
% controller = optimizer(constraints,objective,[],x{1},[u{:}]);
% 
% x = zeros(6,1);
% MPC_x(:,1) = x;
% MPC_u(:,1) = [0;0];
% 
% for i = 1:runtime-1
%     U = controller{x};
%     x = A*x+B*U(:,1);
%    
%    MPC_u(:,i+1) = U(:,1);
%    MPC_x(:,i+1) = x;
%    i
% end
% 
% y = C*MPC_x+D*MPC_u;
% 
% % figure
% % hold on
% % plot(MPC_x(1,:))
% % plot(MPC_x(2,:))
% % plot(MPC_x(3,:))
% % plot(MPC_x(4,:))
% % plot(MPC_x(5,:))
% % plot(MPC_x(6,:))
% % hold off
% 
% figure
% hold on
% plot(y(1,:))
% plot(y(2,:))
% hold off
% 
% figure
% hold on
% plot(MPC_u(1,:))
% plot(MPC_u(2,:))
% hold off

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


