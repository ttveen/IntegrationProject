%% Setup
clear; clc; close all;

%% Load data

addpath('../')

% MPC sin experiment 1
MPCper1 = load('../../Data/MPCper1.mat'); close all;
y_MPCper1 = MPCper1.yhat;
runtime = length(y_MPCper1);
y_MPCper1_r = [35+3*sin((0:1:runtime-1)/100); 30+3*cos((0:1:runtime-1)/150)];   

% PP experiment 2
N = length(y_MPCper1);
PPexp2 = load('../../Data/PPexp2.mat'); close all;
y_pp2_r = PPexp2.r(:,1:N);
y_pp2 = PPexp2.yhat(:,1:N);

% Block Response 5
block = load('../../Data/blockResponse5.mat'); close all;
time = block.time.long(1:N);

%% figures

compSinC = figure('Name','Reference tracking comparison between PP and MPC');
sgtitle('Reference tracking comparison between PP and MPC')

subplot(2,1,1)
title('Temperature Heater 1')
xlabel('Time in s')
ylabel('Temperature in $^{\circ}C$', 'Interpreter', 'Latex')
ylim([0 50]);

hold on
% plot y data and reference
plot(time,y_pp2(1,:),'r+','MarkerIndices',1:20:length(y_pp2(2,:)));
plot(time,y_MPCper1(1,:),'b');
plot(time,y_pp2_r(1,:),'r--');
plot(time,y_MPCper1_r(1,:),'b:');
hold off

legend({'y_{1,PP}','y_{1,MPC}'},'Location','southeast')


subplot(2,1,2)
title('Temperature Heater 2')
xlabel('Time in s')
ylabel('Temperature in $^{\circ}C$', 'Interpreter', 'Latex')
ylim([0 50]);

hold on
% plot y data and reference 
plot(time,y_pp2(2,:),'r+','MarkerIndices',1:20:length(y_pp2(2,:)));
plot(time,y_MPCper1(2,:),'b')
plot(time,y_pp2_r(2,:),'r--');
plot(time,y_MPCper1_r(2,:),'b:');
hold off

legend({'y_{2,PP}','y_{2,MPC}'},'Location','southeast')

% save figure
compSinC = gcf;
compSinC.Renderer = 'painters';
saveas(compSinC, '..\..\Latex\images\controller\compSinC', 'svg');