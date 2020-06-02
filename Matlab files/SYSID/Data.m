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

cd curPath
cd '../..'
cd 'Latex/images/SYSID'
saveas(gcf, 'stepResponse', 'svg');
%% Create an input signal
%Create a block wave
runtime = 200; %200 second runtime
h = 0.50; %sample time
t = 0:h:runtime;
u1 = 0.5*square(t,20)/2+1/2; %Duty cycle of 20%
u2 = -u1
plot(t,u)
%Measure the corresponding output
for i = 1:runtime
    
end