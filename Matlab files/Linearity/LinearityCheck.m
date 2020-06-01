clear all; close all;
%% Check how the system behaves in terms of linearity

% Add some paths
addpath('..\')
addpath('..\')

% include tclab.m for initialization
tclab;

%% First check, input linearity
%f(t,x,u), f(t,x,u) = f(t,x,0.5*u) + f(t,x,0.5*u);

%Ensure that the state is a certain value:
x1 = T1C();
x2 = T2C();

endtime = 100;
h = 0.1;
time1 = 0;
i = 1;

figure(1)

u = 50;

while time1(end) <= endtime
    tic
    %disp("Turn on heaters")
    h1(u);
    h2(u);
    
    t1(i) = T1C();
    t2(i) = T2C();
    
    plot(time1(end),t1(i),'r.')
    hold on;
    plot(time1(end),t2(i),'b.')
    legend('Heater 1', 'Heater 2')
    xlabel('Time in s')
    ylabel('Temperature in /degC',  'Interpreter', 'Latex')
    title('Input at 50%')
    pause(h)
    i = i+1;
    t = toc;
    time1(end+1) = [time1(end) + t];
end
time1(end) = [];
%Turn off heaters
turnOff(a)

%% Wait for the heaters to cool down
t1wait = T1C();
t2wait = T2C();
while (abs(t1wait-x1) > 1 ||abs(t2wait()-x2) > 1)
    t1wait = T1C()
    t2wait = T2C()
    pause(5);
    disp("Waiting")
end
u = 25;
time2 = 0;
i=1;
figure(4)

while time2(end) <= endtime
    tic
    %disp("Turn on heaters")
    h1(u);
    h2(u);
    
    t1_2(i) = T1C();
    t2_2(i) = T2C();
    
    plot(time2(end),t1_2(i),'r.')
    hold on
    plot(time2(end),t2_2(i),'b.')
    legend('Heater 1', 'Heater 2')
    xlabel('Time in s')
    ylabel('Temperature in $^{\circ}C$C',  'Interpreter', 'Latex')
    title('Input at 50%')
    pause(h)
    i = i+1;
    t = toc;
    time2(end+1) = [time2(end) + t];
end
time2(end) =[];
%Turn off heaters
turnOff(a)

%% Make a comparison
t1new = (t1_2-t1wait)*2 + t1wait;
t2new = (t2_2-t2wait)*2 + t2wait;
figure(3)
plot(time2,t1new,'r.')
hold on
plot(time2,t2new,'b.')
legend('Heater 1', 'Heater 2')
xlabel('Time in s')
ylabel('Temperature in /degC',  'Interpreter', 'Latex')
title('Input at 25, doubled%')

figure(4)
plot(time1,t1,'r.', time1, t2, 'b.')
hold on
plot(time2,t1new, 'g.', time2, t2new, 'k.')
legend({'Heater 1 50%', 'Heater 2 50%', 'Heater 1 2*25%', 'Heater 2 2*25%'},'Location', 'northwest')
xlabel('Time in s')
ylabel('Temperature in $^{\circ}C$',  'Interpreter', 'Latex')
title('Comparison: two times half input')
xlabel('Time in s')
ylabel('Temperature in $^{\circ}C$',  'Interpreter', 'Latex')
title('Comparison of input linearity')
cd ..\..
cd 'Latex/Images/Linearity'
saveas(gcf,'Comparison','svg')
%% 2nd Linearity Check
%f(t,x,u), f(t,x,u) = f(t,x,u1) + f(t,x,0.5*u2);
%This means that we compare the heating of the two heaters seperate to
%heating them at the same time

t1wait = T1C();
t2wait = T2C();
while (abs(t1wait-x1) > 1 ||abs(t2wait()-x2) > 1)
    t1wait = T1C()
    t2wait = T2C()
    pause(5);
    disp("Waiting")
end
u = 50;
time3 = 0;
i=1;
figure(5)

while time3(end) <= endtime
    tic
    %disp("Turn on heaters")
    h1(u);
    h2(0);
    
    t1_3(i) = T1C();
    t2_3(i) = T2C();
    
    plot(time3(end),t1_3(i),'r.')
    hold on
    plot(time3(end),t2_3(i),'b.')
    legend('Heater 1', 'Heater 2')
    xlabel('Time in s')
    ylabel('Temperature in $^{\circ}C$C',  'Interpreter', 'Latex')
    title('Input Heater 1 at 50%')
    pause(h)
    i = i+1;
    t = toc;
    time3(end+1) = [time3(end) + t];
end
time3(end) =[];
%Turn off heaters
turnOff(a)

t1wait = T1C();
t2wait = T2C();
while (abs(t1wait-x1) > 1 ||abs(t2wait()-x2) > 1)
    t1wait = T1C()
    t2wait = T2C()
    pause(5);
    disp("Waiting")
end
u = 50;
time4 = 0;
i=1;
figure(6)

while time4(end) <= endtime
    tic
    %disp("Turn on heaters")
    h1(0);
    h2(u);
    
    t1_4(i) = T1C();
    t2_4(i) = T2C();
    
    plot(time4(end),t1_4(i),'r.')
    hold on
    plot(time4(end),t2_4(i),'b.')
    legend('Heater 1', 'Heater 2')
    xlabel('Time in s')
    ylabel('Temperature in $^{\circ}C$C',  'Interpreter', 'Latex')
    title('Input Heater 1 at 50%')
    pause(h)
    i = i+1;
    t = toc;
    time4(end+1) = [time4(end) + t];
end
time4(end) =[];
%Turn off heaters
turnOff(a)
%% Make the 2nd comparison
%Since the index vectors of the two seperate expiriments is not guaranteed
%to be the same length, one measuremnt is interpolated to match the time of
%the measurement with the most number of measurements
timeNew = time3;
if length(time3)>length(time4)
   t1_4 = interp1(time4,t1_4,time3);
   t2_4 = interp1(time4,t2_4,time3);
   timeNew = time3;
elseif length(time3)<length(time4)
   t1_3 = interp1(time3,t1_3,time4);
   t2_3 = interp1(time3,t2_3,time4);
   timeNew = time4;
end

t1new2 = t1_3+t1_4-x1;
t2new2 = t2_3+t2_4-x2;
figure(7)
plot(time1,t1,'r.', time1, t2, 'b.')
hold on
plot(timeNew,t1new2, 'g.', timeNew, t2new2, 'k.')
legend({'Heater 1 50%', 'Heater 2 50%', 'Heater 1 Seperate', 'Heater 2 Seperate'},'Location', 'northwest')
xlabel('Time in s')
ylabel('Temperature in $^{\circ}C$',  'Interpreter', 'Latex')
title('Comparison of input linearity')
xlabel('Time in s')
ylabel('Temperature in $^{\circ}C$',  'Interpreter', 'Latex')
title('Comparison: sum of the two seperate inputs')
cd ..\..
cd 'images/Linearity'
saveas(gcf,'ComparisonSeperate','svg')