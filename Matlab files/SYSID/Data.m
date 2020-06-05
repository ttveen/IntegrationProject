%% Create the experimental data used for system ID
%Add some paths
curPath = pwd;
addpath('../')
% include tclab.m for initialization
clear a;

tclab;
%% Obtain sampling frequency
%Use the step response to get the rise time
runtime = 650; %650 second runtime
figure('Name','Step Response')
fig1a = subplot(2,1,1);
hold on
ylabel('Temperature in $^{\circ}C$',  'Interpreter', 'Latex')
title('Heater 1')
fig1b = subplot(2,1,2);
xlabel('Time in s')
ylabel('Temperature in $^{\circ}C$',  'Interpreter', 'Latex')
title('Heater 2')
hold on

time.step = 0;
for i = 1:runtime
    currenttime = clock;
    tic;
    
    led(0.5);
    h1(50);
    h2(50);

    t1.step(i) = T1C();
    t2.step(i) = T2C();
    

    plot(fig1a,time.step(end),t1.step(i),'r.')
    legend(fig1a,'Heater 1','Location','northwest')
   
    plot(fig1b,time.step(end),t2.step(i),'b.')
    legend(fig1b,'Heater 2','Location','northwest')
    t(i) = toc;
    pause(max(0.01,0.5-t(i)));
    time.step(end+1) = [time.step(end) + etime(clock,currenttime)];
end
sgtitle('Step response at 50%');
hold off

turnOff(a)
led(0)

%cd curPath
cd '../..'
cd 'Latex/images/SYSID'
stepRes = gcf;
stepRes.Renderer = 'painters';
saveas(stepRes, 'stepResponse', 'svg');
%% Create an input signal
%Create a block wave
runtime = 650; %650 second runtime
h = 0.50; %sample time
%t = 0:h:runtim
e;% u1 = 50*square(t,20)/2+1/2; %Duty cycle of 20%
% u2 = -u1
%plot(t,u)
time.block = 0;

figure('Name','Bock wave Response')
fig2a = subplot(2,2,1);
hold on
xlabel('Time in s')
ylabel('Temperature in $^{\circ}C$',  'Interpreter', 'Latex')
title('Heater 1')
fig2b = subplot(2,2,2);
xlabel('Time in s')
ylabel('Temperature in $^{\circ}C$',  'Interpreter', 'Latex')
title('Heater 2')
hold on
fig2c = subplot(2,2,3);
xlabel('Time in s')
ylabel('Input in percentages')
title('Input Heater 1')
hold on
fig2d = subplot(2,2,4);
xlabel('Time in s')
ylabel('Input in percentages')
title('Input Heater 2')
hold on
u1 = zeros(1,runtime);
u2 = zeros(1,runtime);
time.block = zeros(1,runtime);
t = zeros(1,runtime);
t1.block = zeros(1,runtime);
t2.block = zeros(1,runtime);
%Measure the corresponding output
for i = 1:runtime
    currenttime = clock;
    tic
    if time.block < 10
        u1(i) = 0;
        u2(i) = 0;
        led(1);
        h1(50);
        h2(0);
%     elseif time.block < 60
%         u1(i) = 50;
%         u2(i) = 0;
%         led(1);
%         h1(50);
%         h2(0);
    elseif time.block < 120
        u1(i) = 50;
        u2(i) = 0;
        led(1);
        h1(50);
        h2(50)
    elseif time.block < 200
        u1(i) = 10;
        u2(i) = 50;
        led(1);
        h1(10);
        h2(50);
    elseif time.block < 260
        u1(i) = 25;
        u2(i) = 15;
        led(1);
        h1(25);
        h2(15);
    elseif time.block < 350
        u1(i) = 0;
        u2(i) = 0;
        led(1);
        h1(0);
        h2(0);
    elseif time.block < 400
        u1(i) = 25;
        u2(i) = 50;
        led(1);
        h1(25);
        h2(50);
    elseif time.block < runtime
        u1(i) = 50;
        u2(i) = 25;
        led(1);
        h1(50);
        h2(25);
    else
        u1 = 0;
        u2 = 0;
        led(0);
        h1(0);
        h2(0);
    end
    t1.block(i) = T1C();
    t2.block(i) = T2C();
    
    plot(fig2a,time.block(end),t1.block(i),'r.')
%    legend(fig2a,'Heater 1','Location','northwest')
   
    plot(fig2b,time.block(end),t2.block(i),'b.')
%    legend(fig2b,'Heater 2','Location','northwest')
    
    plot(fig2c, time.block(end), u1(i), 'r.')
    
    plot(fig2d, time.block(end), u2(i), 'b.')
    
    t(i) = toc;
    pause(max(0.01,1-t(i)));
    time.block(end+1) = time.block(end) + etime(clock,currenttime);
end

time.block(end) = [];
sgtitle('Response to changing input');
hold off

turnOff(a)
led(0)

%cd 'curPath'
cd '../..'
%cd 'Latex/images/SYSID'
stepBlock = gcf;
stepBlock.Renderer = 'painters';
saveas(stepBlock, 'Latex/images/SYSID/blockResponse2', 'svg');
save('Data/blockResponse2.mat', 'time', 't1', 't2', 'stepBlock', 'u1', 'u2')
%% MUCH DATA
%% Create an input signal
%Create a block wave
runtime = 2*60*60; %2 hour runtime
interval = 100;
intervals = runtime/interval;
time.long = 0;

%Create input
u1 = zeros(1,runtime);
u2 = zeros(1,runtime);
for i = 0:interval
    u1(interval*i+1:interval*(i+1)) = 50*rand(1,1);
    u2(interval*i+1:interval*(i+1)) = 50*rand(1,1);
end
figure('Name','Inputs')
plot(u1); hold on; plot(u2)
%%
figure('Name','Block wave Response, long experiment')
fig3a = subplot(2,1,1);
hold on
xlabel('Time in s')
ylabel('Temperature in ^$\{circ}C$',  'Interpreter', 'Latex')
title('Temperature')
fig3b = subplot(2,1,2);
xlabel('Time in s')
ylabel('Input level',  'Interpreter', 'Latex')
title('Input')
hold on


runtime = 2*60*60; %2 hour runtime
t1.long = zeros(1,runtime);
t2.long = zeros(1,runtime);
t = zeros(1,runtime);
time.long = zeros(1,runtime);
%Measure the corresponding output
led(1);
try
for i = 1:runtime
    currenttime = clock;
    tic
    h1(u1(i));
    h2(u2(i));
    t1.long(i) = T1C();
    t2.long(i) = T2C();
    if mod(i,100) == 0
       disp(t(i))
       disp(time.long(i))
       disp(getTemp())
    end
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

time.long(end) = [];
plot(fig3a,time.long,t1.long,'r.')
plot(fig3a,time.long,t2.long,'b.')
legend(fig3a,{'Heater 1','Heater 2'},'Location','northwest')

plot(fig3b, time.long, u1(1:length(time.long)), 'r.')
plot(fig3b, time.long, u2(1:length(time.long)), 'b.')
legend(fig3b,{'Input 1', 'Input 2'},'Location','northwest')


sgtitle('Response to changing input');
hold off



%cd 'curPath'
cd '../..'
%cd 'Latex/images/SYSID'
stepBlock = gcf;
stepBlock.Renderer = 'painters';
saveas(stepBlock, 'Latex/images/SYSID/blockResponse5', 'svg');
save('Data/blockResponse5.mat', 'time', 't1', 't2', 'stepBlock', 'u1', 'u2')