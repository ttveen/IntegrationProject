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


%% Calculate errors
N = size(y_exp,2);

% root mean squared error (RMSE)
RMSE_moesp_y1 = sqrt((1/N)*(sum((y_exp(1,:) - y_moesp(1,:)).^2)))
RMSE_moesp_y2 = sqrt((1/N)*(sum((y_exp(2,:) - y_moesp(2,:)).^2)))

n = 400;
RMSE_nsid_y1 = sqrt((1/N)*(sum((y_exp(1,1:n) - y_nsid(1,1:n)).^2)))
RMSE_nsid_y2 = sqrt((1/N)*(sum((y_exp(2,1:n) - y_nsid(2,1:n)).^2)))


