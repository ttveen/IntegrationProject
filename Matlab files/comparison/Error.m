%% Setup
clear; clc; close all;

%% Load Data Sets
addpath('../')

% Block Response 5
load('../../Data/blockResponse5.mat'); close all;
y_exp = [t1.long; t2.long];
u_exp = [u1; u2];
time = time.long;

% SID
load('../../Data/SID.mat'); close all;
y_moesp = SID.ys; 

% NSID
load('../../Data/NSID.mat'); close all;
y_nsid = NSID.ys;

% Least squares
LS = load('../../Data/LeastSquares.mat'); close all;
LS.TC1 = interp1(LS.t,LS.TC1,time);
LS.TC2 = interp1(LS.t,LS.TC2,time);
y_ls = [LS.TC1; LS.TC2]; 

% PP experiment 1
N = 1000;
PPexp1 = load('../../Data/PPexp1.mat'); close all;
y_pp1 = PPexp1.yhat(:,1:N);
y_pp1_r = 30.*ones(2,length(y_pp1));

% PP experiment 2
PPexp2 = load('../../Data/PPexp2.mat'); close all;
y_pp2_r = PPexp2.r(:,1:N);
y_pp2 = PPexp2.yhat(:,1:N);

% Kalman filter experiment 1
KF1 = load('../../Data/kalmanTest1.mat'); close all;
y_kf1_exp = [KF1.t1.kalman; KF1.t2.kalman];
y_kf1 = KF1.yhat;

% Kalman filter experiment 2
KF2 = load('../../Data/kalmanTest2.mat'); close all;
y_kf2_exp = [KF2.t1.kalman; KF2.t2.kalman];
y_kf2 = KF2.yhat;

% MPC step ref
MPCstep = load('../../Data/MPCstep.mat'); close all;
y_MPCstep = MPCstep.yhat; 
y_MPCstep_r = 30.*ones(2,length(y_MPCstep));

% MPC sin experiment 1
MPCper1 = load('../../Data/MPCper1.mat'); close all;
y_MPCper1 = MPCper1.yhat;
runtime = length(y_MPCper1);
y_MPCper1_r = [35+3*sin((0:1:runtime-1)/100); 30+3*cos((0:1:runtime-1)/150)];   

% MPC sin experiment 2
MPCper2 = load('../../Data/MPCper2.mat'); close all;
y_MPCper2 = MPCper2.yhat;
runtime = length(y_MPCper2);
y_MPCper2_r = [35+10*sin((0:1:runtime-1)/50); 30+5*cos((0:1:runtime-1)/75)];

% LQR
LQR = load('../../Data/LQRexp1.mat'); close all;
N = 1243;
y_lqr = LQR.yhat(:,1:N);
y_lqr_r = 30*ones(size(y_lqr));
%% Calculate errors

% root mean squared error (RMSE)
% moesp (output)
RMSE_moesp_y1 = rmse(y_exp(1,:),y_moesp(1,:));
RMSE_moesp_y2 = rmse(y_exp(2,:),y_moesp(2,:));
RMSE_moesp = [RMSE_moesp_y1; RMSE_moesp_y2]

% n4sid (output)
RMSE_nsid_y1 = rmse(y_exp(1,:),y_nsid(1,:));
RMSE_nsid_y2 = rmse(y_exp(2,:),y_nsid(2,:));

n = 400;
RMSE_nsid_y1n = rmse(y_exp(1,n:end),y_nsid(1,n:end));
RMSE_nsid_y2n = rmse(y_exp(2,n:end),y_nsid(2,n:end));
RMSE_nsid = [RMSE_nsid_y1 RMSE_nsid_y1n; RMSE_nsid_y2 RMSE_nsid_y2n]

% least squares (output)
RMSE_ls_y1 = rmse(y_exp(1,:),y_ls(1,:));
RMSE_ls_y2 = rmse(y_exp(2,:),y_ls(2,:));
RMSE_ls = [RMSE_ls_y1; RMSE_ls_y2]

% pole placement experiment 1 (ref)
RMSE_pp1_y1 = rmse(y_pp1_r(1,:),y_pp1(1,:));
RMSE_pp1_y2 = rmse(y_pp1_r(2,:),y_pp1(2,:));
RMSE_pp1 = [RMSE_pp1_y1; RMSE_pp1_y2]

% pole placement experiment 2 (ref)
RMSE_pp2_y1 = rmse(y_pp2_r(1,:),y_pp2(1,:));
RMSE_pp2_y2 = rmse(y_pp2_r(2,:),y_pp2(2,:));
RMSE_pp2 = [RMSE_pp2_y1; RMSE_pp2_y2]

% kalman filter experiment 1 (output) 
RMSE_kf1_y1 = rmse(y_kf1_exp(1,:),y_kf1(1,:));
RMSE_kf1_y2 = rmse(y_kf1_exp(2,:),y_kf1(2,:));

n = 40;
RMSE_kf1_y1n = rmse(y_kf1_exp(1,n:end),y_kf1(1,n:end));
RMSE_kf1_y2n = rmse(y_kf1_exp(2,n:end),y_kf1(2,n:end));
RMSE_kf1 = [RMSE_kf1_y1 RMSE_kf1_y1n; RMSE_kf1_y2 RMSE_kf1_y2n]

% kalman filter experiment 2 (output) 
RMSE_kf2_y1 = rmse(y_kf2_exp(1,:),y_kf2(1,:));
RMSE_kf2_y2 = rmse(y_kf2_exp(2,:),y_kf2(2,:));

n = 40;
RMSE_kf2_y1n = rmse(y_kf2_exp(1,n:end),y_kf2(1,n:end));
RMSE_kf2_y2n = rmse(y_kf2_exp(2,n:end),y_kf2(2,n:end));
RMSE_kf2 = [RMSE_kf2_y1 RMSE_kf2_y1n; RMSE_kf2_y2 RMSE_kf2_y2n]

% MPC step (ref)
RMSE_MPCstep_y1 = rmse(y_MPCstep_r(1,:),y_MPCstep(1,:));
RMSE_MPCstep_y2 = rmse(y_MPCstep_r(2,:),y_MPCstep(2,:));
RMSE_MPCstep = [RMSE_MPCstep_y1; RMSE_MPCstep_y2]

% MPC sin experiment 1
RMSE_MPCper1_y1 = rmse(y_MPCper1_r(1,:),y_MPCper1(1,:));
RMSE_MPCper1_y2 = rmse(y_MPCper1_r(2,:),y_MPCper1(2,:));
RMSE_MPCper1 = [RMSE_MPCper1_y1; RMSE_MPCper1_y2]

% MPC sin experiment 2
RMSE_MPCper2_y1 = rmse(y_MPCper2_r(1,:),y_MPCper2(1,:));
RMSE_MPCper2_y2 = rmse(y_MPCper2_r(2,:),y_MPCper2(2,:));
RMSE_MPCper2 = [RMSE_MPCper2_y1; RMSE_MPCper2_y2]

% LQR
RMSE_lqr_y1 = rmse(y_lqr_r(1,:),y_lqr(1,:));
RMSE_lqr_y2 = rmse(y_lqr_r(2,:),y_lqr(2,:));
RMSE_lqr = [RMSE_lqr_y1; RMSE_lqr_y2]

%% Keep those pesky figures out of here
close all;

%% functions
function rms = rmse(y1,y2)
    N = length(y1);
    rms = sqrt((1/N)*(sum((y1 - y2).^2)));
end