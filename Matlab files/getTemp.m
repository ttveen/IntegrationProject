function [T1, T2] = getTemp()
% getTemp returns the temperature of the two heaters

% include tclab.m for initialization
tclab;

% get temperature and display it in the command window.
T1 = T1C();
T2 = T2C();

disp(['Temperature 1: ' num2str(T1) ' degC'])
disp(['Temperature 2: ' num2str(T2) ' degC'])
end

