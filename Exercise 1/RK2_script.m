t_i = 0; % set initial value of t_0
q_i = 500e-9; % set q_initial condition q at t_0
h = 0.0003; % set t step-size
t_f = 1; % stop here
R = 1000; % resistance
C = 100e-9; % capacitance


%******************Step_Signal 2.5V******************
func = @(t, q) 2.5/R - 1/(R*C)*q; %function handle: 2 variables
%****************************************************


% %******************Impulsive_Signal 2.5V******************
% tau = 100;
% func = @(t, q) 2.5/R*exp(-t^2/tau) - 1/(R*C)*q; %function handle: 2 variables
% %*********************************************************


% %******************Decay_Signal 2.5V******************
% tau = 100;
% func = @(t, q) 2.5/R*exp(-t/tau) - 1/(R*C)*q; %function handle: 2 variables
% %*****************************************************


% %******************Sine_Wave_Signal******************
% period = 10e-6;
% period = 100e-6;
% period = 500e-6;
% period = 1000e-6;
% func = @(t, q) 5/R*sin(2*pi/period*t) - 1/(R*C)*q; %function handle: 2 variables
% %*****************************************************


% %******************Squre_Wave_Signal******************
% period = 10e-6;
% period = 100e-6;
% period = 500e-6;
% period = 1000e-6;
% func = @(t, q) 5/R*square(2*pi/period*t) - 1/(R*C)*q; %function handle: 2 variables
% %*****************************************************


% %******************Sawtooth_Wave_Signal******************
% period = 10e-6;
% period = 100e-6;
% period = 500e-6;
% period = 1000e-6;
% func = @(t, q) 5/R*sawtooth(2*pi/period*t) - 1/(R*C)*q; %function handle: 2 variables
% %********************************************************


[xout, yout] = RK2(func, t_i, q_i, t_f, h); % call to euler.m
plot(xout, yout, 'b') %plot now values of x,y