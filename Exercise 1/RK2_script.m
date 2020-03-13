t_i = 0; % set initial value of t_0
q_i = 500e-9; % set q_initial condition q at t_0
h = 0.00001; % set t step-size
t_f = 0.01; % stop here
R = 1000; % resistance
C = 100e-9; % capacitance


%******************Step_Signal 2.5V******************
% func = @(t, q) 2.5/R*heaviside(t) - 1/(R*C)*q; %function handle: 2 variables
%
% [tout, qout] = RK2(func, t_i, q_i, t_f, h); % call to euler.m
% plot(tout, qout, 'b') %plot now values of x,y
% 
% figure; % use this if you want the exact solution and the numerical on different figures
% %hold on; % use this if you want the exact solution and the numerical on the same figure
% 
% qexact = 2.5 * C * ( 1 + exp(-tout/(R*C)) ); %calculate exact solution
% plot(tout, qexact, 'r') %plot now values of x,y
%****************************************************


% %******************Impulsive_Signal 2.5V******************
% tau = 100;
% func = @(t, q) 2.5/R*exp(-t^2/tau) - 1/(R*C)*q; %function handle: 2 variables
% 
% [tout, qout] = RK2(func, t_i, q_i, t_f, h); % call to euler.m
% plot(tout, qout, 'b') %plot now values of x,y
% 
% There cannot be an exact solution for this Vin
% %*********************************************************


% %******************Decay_Signal 2.5V******************
% tau = 100;
% func = @(t, q) 2.5/R*exp(-t/tau) - 1/(R*C)*q; %function handle: 2 variables
% 
% [tout, qout] = RK2(func, t_i, q_i, t_f, h); % call to euler.m
% plot(tout, qout, 'b') %plot now values of x,y
% 
% figure; % use this if you want the exact solution and the numerical on different figures
% % %hold on; % use this if you want the exact solution and the numerical on the same figure
% 
% qexact = 2.5 * tout * exp(-tout/(R*C)) / R + 5 * C * exp(-tout/(R*C)); %calculate exact solution
% plot(tout, qexact, 'r') %plot now values of x,y
% %*****************************************************


% %******************Sine_Wave_Signal******************
% period = 10e-6;
% period = 100e-6;
% period = 500e-6;
% period = 1000e-6;
% func = @(t, q) 5/R*sin(2*pi/period*t) - 1/(R*C)*q; %function handle: 2 variables
% 
% [tout, qout] = RK2(func, t_i, q_i, t_f, h); % call to euler.m
% plot(tout, qout, 'b') %plot now values of x,y
% 
% figure; % use this if you want the exact solution and the numerical on different figures
% %hold on; % use this if you want the exact solution and the numerical on the same figure
% 
% qexact = 2.5*C*period/(period^2 + (2*pi*R*C)^2)*(period*sin(2*pi/period*tout) - 2*pi*C*R*cos(2*pi/period*tout)) + (q_i + (2.5*C^2*period*2*pi*R)/(period^2 + (2*pi*C*R)^2))*exp(-tout/(R*C)); %calculate exact solution
% plot(tout, qexact, 'r') %plot now values of x,y
% %*****************************************************


% %******************Squre_Wave_Signal******************
% period = 10e-6;
% period = 100e-6;
% period = 500e-6;
% period = 1000e-6;
% func = @(t, q) 5/R*square(2*pi/period*t) - 1/(R*C)*q; %function handle: 2 variables
% 
% [tout, qout] = RK2(func, t_i, q_i, t_f, h); % call to euler.m
% plot(tout, qout, 'b') %plot now values of x,y
% 
% There cannot be an exact solution for this Vin
% %*****************************************************


% %******************Sawtooth_Wave_Signal******************
% period = 10e-6;
% period = 100e-6;
% period = 500e-6;
% period = 1000e-6;
% func = @(t, q) 5/R*sawtooth(2*pi/period*t) - 1/(R*C)*q; %function handle: 2 variables
% 
% [tout, qout] = RK2(func, t_i, q_i, t_f, h); % call to euler.m
% plot(tout, qout, 'b') %plot now values of x,y
% 
% There cannot be an exact solution for this Vin
% %********************************************************