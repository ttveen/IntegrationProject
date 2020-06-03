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
ylabel('Temperature in /degC',  'Interpreter', 'Latex')
title('Heater 1')
fig1b = subplot(2,1,2);
xlabel('Time in s')
ylabel('Temperature in /degC',  'Interpreter', 'Latex')
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
t = 0:h:runtime;
% u1 = 50*square(t,20)/2+1/2; %Duty cycle of 20%
% u2 = -u1
%plot(t,u)
time.block = 0;

figure('Name','Bock wave Response')
fig2a = subplot(2,2,1);
hold on
xlabel('Time in s')
ylabel('Temperature in /degC',  'Interpreter', 'Latex')
title('Heater 1')
fig2b = subplot(2,2,2);
xlabel('Time in s')
ylabel('Temperature in /degC',  'Interpreter', 'Latex')
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

%Measure the corresponding output
for i = 1:runtime
    currenttime = clock;
    tic
    if time.block < 50
        u1(i) = 50;
        u2(i) = 0;
        led(1);
        h1(50);
        h2(0);
    elseif time.block < 120
        u1(i) = 25;
        u2(i) = 50;
        led(1);
        h1(25);
        h2(50)
    elseif time.block < 200
        u1(i) = 10;
        u2(i) = 40;
        led(1);
        h1(10);
        h2(40);
    elseif time.block < 260
        u1(i) = 40;
        u2(i) = 40;
        led(1);
        h1(40);
        h2(40);
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
    pause(max(0.01,0.5-t(i)));
    time.block(end+1) = [time.block(end) + etime(clock,currenttime)];
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
saveas(stepBlock, 'Latex/images/SYSID/blockResponse', 'svg');
save('Data/blockResponse.mat', 'time', 't1', 't2', 'stepBlock', 'u1', 'u2')