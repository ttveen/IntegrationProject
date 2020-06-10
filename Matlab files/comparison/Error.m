%% Setup
clear; clc; close all;

%% Load Data Sets
addpath('../')

% Block Response 5
load('../../Data/blockResponse5.mat');
y_exp = [t1.long; t2.long];
u_exp = [u1; u2];
time = time.long;

% SID
load('../../Data/SID.mat');
y_moesp = SID.ys;

% NSID
load('../../Data/NSID.mat');
y_nsid = NSID.ys;

% % PP experiment 2
PPexp2 = load('../../Data/PPexp2.mat');
y_exp_pp2 = [PPexp2.t1.PP; PPexp2.t2.PP];
y_pp2 = PPexp2.yhat;

% Kalman filter
KF = load('../../Data/kalmanTest2.mat');
y_exp_kf = [KF.t1.kalman; KF.t2.kalman];
y_kf = KF.yhat;

%% Calculate errors
N = size(y_exp,2);

% root mean squared error (RMSE)
RMSE_moesp_y1 = sqrt((1/N)*(sum((y_exp(1,:) - y_moesp(1,:)).^2)));
RMSE_moesp_y2 = sqrt((1/N)*(sum((y_exp(2,:) - y_moesp(2,:)).^2)));
RMSE_moesp = [RMSE_moesp_y1; RMSE_moesp_y2]

n = 400;
RMSE_nsid_y1 = rmse(y_exp(1,1:n),y_nsid(1,1:n));
RMSE_nsid_y2 = rmse(y_exp(2,1:n),y_nsid(2,1:n))
RMSE_nsid = [RMSE_nsid_y1; RMSE_nsid_y2]

RMSE_pp2_y1 = rmse(y_exp_pp2(1,:),y_pp2(1,:));
RMSE_pp2_y2 = rmse(y_exp_pp2(2,:),y_pp2(2,:));
RMSE_pp2 = [RMSE_pp2_y1; RMSE_pp2_y2]

RMSE_kf_y1 = sqrt((1/N)*(sum((y_exp_pp2(1,:) - y_pp2(1,:)).^2)));

function rms = rmse(y1,y2)
    N = length(y1);
    rms = sqrt((1/N)*(sum((y1 - y2).^2)));
end