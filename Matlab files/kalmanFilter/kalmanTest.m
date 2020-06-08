%% Kalman filter test
%%
load('../../Data/N4SID');
tclab;
%%
figure('Name','Kalman Filter Test')
fig1a = subplot(2,1,1);
hold on
xlabel('Time in s')
ylabel('Temperature in ^$\{circ}C$',  'Interpreter', 'Latex')
title('Temperature')
fig1b = subplot(2,1,2);
xlabel('Time in s')
ylabel('Input level',  'Interpreter', 'Latex')
title('Input')
hold on


runtime = 200; %200 second runtime
t1.kalman = zeros(1,runtime);
t2.kalman = zeros(1,runtime);
t = zeros(1,runtime);
time.kalman = zeros(1,runtime);
xhat = zeros(size(A,1),runtime)
yhat = zeros(zeros(2,runtime));

%Create input
interval = 20; %20 different inputs
intervals = runtime/interval;
u1 = zeros(1,runtime);
u2 = zeros(1,runtime);
for i = 0:interval
    u1(interval*i+1:interval*(i+1)) = 50*rand(1,1);
    u2(interval*i+1:interval*(i+1)) = 50*rand(1,1);
end

%Measure the corresponding output
led(1);
try
for i = 1:runtime
    currenttime = clock;
    tic
    h1(u1(i));
    h2(u2(i));
    t1.kalman(i) = T1C();
    t2.kalman(i) = T2C();
    
    %Kalman filter
    xhat(i+1) = (A-K*C)*xhat(i) + B*[u1(i); u2(i)] + K*[t1.kalman(i); t2.kalman(i)];
    yhat(i) = C*xhat(i) + D*[u1(i); u2(i)];
    
    plot(fig1a,time.block(end),t1.block(i),'r.')
    hold on
    plot(fig1a,time.block(end),t2.block(i),'b.')
    %legend(fig2a,'Heater 1','Location','northwest')
    
    plot(fig1b, time.block(end), u1(i), 'r.')
    hold on
    plot(fig2b, time.block(end), u2(i), 'b.')
    %legend(fig2b,'Heater 2','Location','northwest')
    
    t(i) = toc;
    pause(max(0.01,1-t(i)));
    time.long(i+1) = time.long(i) + etime(clock,currenttime);
end
catch
    disp('error')
   turnOff(a)
   led(0)
end
turnOff(a)
led(0)