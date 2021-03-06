%% Kalman filter test
clear all
close all
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
runtime = 2000; %200 second runtime
% axis([fig1a fig1b], [0 runtime 0 50]);
t = zeros(1,runtime);
time.PP = zeros(1,runtime);

% init Kalman filter
xhat = zeros(size(A,1),runtime);
yhat = zeros(2,runtime);

% init measurements
t1.PP = zeros(1,runtime);
t2.PP = zeros(1,runtime);

% init input
u = zeros(2,runtime);

% init reference
r = [35 + 3*sin((0:1:runtime-1)/100); 30 + 3*cos((0:1:runtime-1)/100)];

%% Determine gain
% p = [-0.9 -0.8 -0.7 0.7 0.8 0.9]; % stable discrete poles
% p = [-0.3 -0.2 -0.1 0.1 0.2 0.3];
p = linspace(0.987,0.999,6);


Kx = place(A,B,p); % pole placement
%Kr = inv(C/(eye(size(A))-A+B*Kx)*B); % reference gain
%For LQR control
%Kx = [-22.3112995358586,-2.34395829803522,1.44421673670272,0.463052019339851,-0.586068184551735,0.0519736559802679;-16.2696117036986,2.56721439335664,1.24800644905354,1.20426827073047,-0.460281038324108,0.0545396527165367];
Kr = inv((C-D*Kx)*inv(eye(size(A))-A+B*Kx)*B+D);
% %% simulate
x = xhat;

x(:,1) = (C-D*Kx)\([23;23]-D*Kr*r(:,1)); 
for i = 1:runtime-1
    % generated input
    u(:,i) = -Kx*x(:,i) + Kr*r(:,i);
    
    % system
    x(:,i+1) = A*x(:,i)+B*u(:,i);
        
end
y = C*x+D*u;
close all
figure('Name','Pole Placement Simulation')
fig2a = subplot(2,1,1);
plot(y(1,:))
hold on
plot(y(2,:))
title('Simulation output')
legend('y1','y2')%,'r1','r2')
fig2b = subplot(2,1,2);
plot(u(1,:))
hold on
plot(u(2,:))
title('Controller output / system input')
legend('u1','u2')
sgtitle('Pole Placement step reference tracking simulation')
% PPstepref = gcf;
% saveas(PPstepref, '../../Latex/images/controller/PPsetref1', 'svg');


%% Create figure
figure('Name','PPexp')
fig1a = subplot(2,1,1);
aniT1 = animatedline(fig1a,'Color','r','Marker','.','Linestyle','none');
aniT2 = animatedline(fig1a,'Color','b','Marker','.','Linestyle','none');
aniK1 = animatedline(fig1a,'Color','r','Marker','+','Linestyle','none');
aniK2 = animatedline(fig1a,'Color','b','Marker','+','Linestyle','none');
hold on
xlabel('Time in s')
ylabel({'Temperature in $^{\circ}C$'},  'Interpreter', 'Latex')
title('Temperature')
axis([0 runtime 0 35])
fig1b = subplot(2,1,2)
aniu1 = animatedline(fig1b,'Color','r','Marker','.','Linestyle','none');
aniu2 = animatedline(fig1b,'Color','b','Marker','.','Linestyle','none');
xlabel('Time in s')
ylabel('Input level',  'Interpreter', 'Latex')
title('Input')
axis([0 runtime -5 50])
hold on
plot(fig1a,(0:1:runtime-1),r(1,:),'r--')
plot(fig1a,(0:1:runtime-1),r(2,:),'b--')
%Measure the corresponding output
led(1);
%% Run Experiment
for i = 1:runtime - 1
    currenttime = clock;
    tic
    
    % generated input
    u(:,i) = -Kx*xhat(:,i) + Kr*r(:,i);
    % implement input
    h1(u(1,i));
    h2(u(2,i));
    
    % measure output
    t1.PP(i) = T1C();
    t2.PP(i) = T2C();
    
    % Kalman filter
    xhat(:,i+1) = (A-K*C)*xhat(:,i) + B*[u(1,i); u(2,i)] + K*[t1.PP(i); t2.PP(i)];
    yhat(:,i) = C*xhat(:,i) + D*[u(1,i); u(2,i)];
    
    
    % plot results
    addpoints(aniT1,time.PP(i),t1.PP(i))
    %plot(fig1a,time.PP(i),t1.PP(i),'r.')
    hold on
    addpoints(aniT2,time.PP(i),t2.PP(i))
    addpoints(aniK1,time.PP(i),yhat(1,i))
    addpoints(aniK2,time.PP(i),yhat(2,i))
    %plot(fig1a,time.PP(i),t2.PP(i),'b.')
    %plot(fig1a,time.PP(i), yhat(1,i), 'r+')
    %plot(fig1a,time.PP(i), yhat(2,i), 'b+')
    %legend(fig2a,'Heater 1','Location','northwest')
    
    addpoints(aniu1, time.PP(i), u(1,i))
    hold on
    addpoints(aniu2, time.PP(i), u(2,i))
    %legend(fig2b,'Heater 2','Location','northwest')
    
    % time stuff
    t(i) = toc;
    pause(max(0.01,1-t(i)));
    time.PP(i+1) = time.PP(i) + etime(clock,currenttime);
end

turnOff(a)
led(0)
time.PP(end) = [];

legend(fig1a,{'$T_{C1,reference}$','$T_{C1,reference}$','$T_{C1,measured}$','$T_{C2,measured}$','$T_{C1,estimate}$','$T_{C2,estimate}$'},'Interpreter','latex','Location','southeast','FontSize', 20)
legend(fig1b,{'$u_1$', '$u_2$'}, 'Interpreter','latex','FontSize', 20)
PPexp = gcf;
set(PPexp, 'position', get(0, 'ScreenSize'))
PPexp.Renderer = 'painters';
saveas(PPexp, '../../Latex/images/Controller/PPexp2', 'svg');
save('../../Data/PPexp2.mat', 'time', 't1', 't2', 'PPexp', 'u', 'xhat', 'yhat','r')
close gcf