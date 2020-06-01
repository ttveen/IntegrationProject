clear all; close all;
%% Check how the system behaves in terms of linearity

% include tclab.m for initialization
tclab;

%% First check, input linearity
% f(t,x,u), f(t,x,u) = f(t,x,0.5*u) + f(t,x,0.5*u);

%Ensure that the state is a certain value:
x1 = T1C();
x2 = T2C();

endtime = 100;
h = 0.1;
time1 = 0;
i = 1;

figure(1)
% plot(time,x1,'r.')
% hold on;
% plot(time,x2,'b.')

t = 0:h:endtime;

u = 50;

%t1 = zeros(1,endtime/h);
%t2 = t1;
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
while (abs(t1wait-x1) > 2 ||abs(t2wait()-x2) > 2)
    t1wait = T1C()
    t2wait = T2C()
    pause(5);
    disp("Waiting")
end
u = 25;
time2 = 0;
i=1;
figure(2)
% plot(time,t1wait,'r.')
% hold on;
% plot(time,t2wait,'b.')

xlabel('Time in s')
ylabel('Temperature in $^{\circ}C$C',  'Interpreter', 'Latex')
title('Input at 25%')

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

figure(4)
plot(time1,t1,'r.', time1, t2, 'b.')
hold on
plot(time2,t1new, 'g.', time2, t2new, 'k.')
legend({'Heater 1 50%', 'Heater 2 50%', 'Heater 3 2*25%', 'Heater 3 2*25%'},'Location', 'northwest')
xlabel('Time in s')
ylabel('Temperature in $^{\circ}C$',  'Interpreter', 'Latex')
title('Comparison of input linearity')
xlabel('Time in s')
ylabel('Temperature in $^{\circ}C$',  'Interpreter', 'Latex')
title('Comparison of input linearity')
saveas(gcf,