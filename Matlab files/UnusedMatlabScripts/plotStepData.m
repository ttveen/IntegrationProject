% Temporary script to plot the step responses with the saved data on drive
% The data has to be manually loaded in

figure('Name','Step Response')
sgtitle('Step response at 50%');

fig1a = subplot(2,1,1);
title('Heater 1')
ylabel('Temperature in /degC',  'Interpreter', 'Latex')

plot(fig1a,time.step,t1.step,'r.')
legend(fig1a,'Heater 1','Location','northwest')
% yline(67.79,'b') % SS
% yline(29.96,'b') % 10%
% yline(63.58,'b')% 90%

fig1b = subplot(2,1,2);
title('Heater 2')
xlabel('Time in s')
ylabel('Temperature in /degC',  'Interpreter', 'Latex')

plot(fig1b,time.step,t2.step,'b.')
legend(fig1b,'Heater 2','Location','northwest')
% yline(57.04,'b') % SS
% yline(28.89,'b') % 10%
% yline(53.91,'b')% 90%

