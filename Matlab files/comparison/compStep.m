%% Setup
clear; clc; close all;

%% Load data

addpath('../')

% MPC step ref
MPCstep = load('../../Data/MPCstep.mat'); close all;
y_MPCstep = MPCstep.yhat; 

% PP experiment 1
N = length(y_MPCstep);
PPexp1 = load('../../Data/PPexp1.mat'); close all;
y_pp1 = PPexp1.yhat(:,1:N);

% Block Response 5
load('../../Data/blockResponse5.mat'); close all;
time = time.long(1:N);

r = 30.*ones(2,N);

%% Step info
Spp_y1 = stepinfo(y_pp1(1,:),time,30);
Spp_y2 = stepinfo(y_pp1(2,:),time,30);

Smpc_y1 = stepinfo(y_MPCstep(1,:),time,30);
Smpc_y2 = stepinfo(y_MPCstep(2,:),time,30);

%% figures

compStepC = figure('Name','Comparison PP and MPC (Step)')
sgtitle('Step response comparison between PP and MPC')

subplot(2,1,1)
title('Temperature Heater 1')
xlabel('Time in s')
ylabel('Temperature in $^{\circ}C$', 'Interpreter', 'Latex')
ylim([0 35]);

hold on
% plot y data and reference
plot(time,y_pp1(1,:),'r+','MarkerIndices',1:20:length(y_pp1(1,:)));
plot(time,y_MPCstep(1,:),'b');
plot(time,r(1,:),'k:');

% plot stepinfo
hold off

legend({'y_{1,PP}','y_{1,MPC}'},'Location','southeast')


subplot(2,1,2)
title('Temperature Heater 2')
xlabel('Time in s')
ylabel('Temperature in $^{\circ}C$', 'Interpreter', 'Latex')
ylim([0 35]);

hold on
% plot y data and reference 
plot(time,y_pp1(2,:),'r+','MarkerIndices',1:20:length(y_pp1(2,:)));
plot(time,y_MPCstep(2,:),'b')
plot(time,r(2,:),'k:');
hold off

legend({'y_{2,PP}','y_{2,MPC}'},'Location','southeast')

% save figure
compStepC = gcf;
compStepC.Renderer = 'painters';
saveas(compStepC, '..\..\Latex\images\controller\compStepC', 'svg');