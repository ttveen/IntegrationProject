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
N = 30; % seconds

% init reference
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