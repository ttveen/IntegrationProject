%% Kalman filter test
%%
addpath('../')
load('../../Data/NSID');
tclab;
A = NSID.At;
B = NSID.Bt;
C = NSID.Ct;
D = NSID.Dt;
K = NSID.Ks;
%%
figure('Name','Kalman Filter Test')
fig1a = subplot(2,1,1);
hold on
xlabel('Time in s')
ylabel({'Temperature in $^{\circ}C$'},  'Interpreter', 'Latex')
title('Temperature')
fig1b = subplot(2,1,2);
xlabel('Time in s')
ylabel('Input level',  'Interpreter', 'Latex')
title('Input')
hold on


runtime = 200; %200 second runtime
axis([fig1a fig1b], [0 runtime 0 50]);
t1.kalman = zeros(1,runtime);
t2.kalman = zeros(1,runtime);
t = zeros(1,runtime);
time.kalman = zeros(1,runtime);
xhat = zeros(size(A,1),runtime);
yhat = zeros(2,runtime);

%Create input
interval = 50; %50 second interval
intervals = runtime/interval;
u1 = zeros(1,runtime);
u2 = zeros(1,runtime);
for i = 0:interval
    u1(interval*i+1:interval*(i+1)) = 50*rand(1,1);
    u2(interval*i+1:interval*(i+1)) = 50*rand(1,1);
end

%Measure the corresponding output
led(1);
%try
for i = 1:runtime
    currenttime = clock;
    tic
    h1(u1(i));
    h2(u2(i));
    t1.kalman(i) = T1C();
    t2.kalman(i) = T2C();
    
    %Kalman filter
    xhat(:,i+1) = (A-K*C)*xhat(:,i) + B*[u1(i); u2(i)] + K*[t1.kalman(i); t2.kalman(i)];
    yhat(:,i) = C*xhat(:,i) + D*[u1(i); u2(i)];
    
    plot(fig1a,time.kalman(i),t1.kalman(i),'r.')
    hold on
    plot(fig1a,time.kalman(i),t2.kalman(i),'b.')
    plot(fig1a,time.kalman(i), yhat(1,i), 'r+')
    plot(fig1a,time.kalman(i), yhat(2,i), 'b+')
    %legend(fig2a,'Heater 1','Location','northwest')
    
    plot(fig1b, time.kalman(i), u1(i), 'r.')
    hold on
    plot(fig1b, time.kalman(i), u2(i), 'b.')
    %legend(fig2b,'Heater 2','Location','northwest')
    
    t(i) = toc;
    pause(max(0.01,1-t(i)));
    time.kalman(i+1) = time.kalman(i) + etime(clock,currenttime);
end
% catch
%     disp('error')
%    turnOff(a)
%    led(0)
% end
turnOff(a)
led(0)
time.kalman(end) = [];

kalmanTestFig = gcf;
stepBlock.Renderer = 'painters';
saveas(kalmanTestFig, '../../Latex/images/kalmanTest/kalmanTest2', 'svg');
save('../../Data/kalmanTest2.mat', 'time', 't1', 't2', 'kalmanTestFig', 'u1', 'u2', 'xhat', 'yhat')