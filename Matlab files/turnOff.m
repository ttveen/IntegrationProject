function turnOff()
% turnOff turns of the two heaters

% include tclab.m for initialization
tclab;

% turn off the heaters by setting both inputs to zero
h1(0); % heater 1
h2(0); % heater 2

disp('The heaters have been turned off');
end

