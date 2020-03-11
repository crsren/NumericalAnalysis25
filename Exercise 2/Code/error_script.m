R = 1000; %resistance.
C = 100/1000000000; %capacitance.
qc0 = 500/1000000000; %charge at t0.
t0 = 0; %value of time when time is zero (obviously).
h = 0.00001; %step size.
tf = 0.01; %interval [0, 0.01].

T = 100/1000000; %period of T (1/T = frequency). 2pi*frequency = rotationshstighet (angular frequency?).
freq = 1/T; % frequency 
w = pi * freq; %angular frequency
V = 5.0; %amplitude (voltage).


%****************Heun's Method****************
%      a = 0.5;
%*********************************************

%***************Midpoint Method***************
%      a = 0.0;
%*********************************************

%*****************Random RK2******************
%      a = 0.3;
%*********************************************



Vi = @(t) V * cos(w*t); %input cosine wave with given values.

qcFunc = @(t, qc) ( Vi(t)/R ) - ( 1/(R*C) ) * qc; %we get the d/dt of qc (derivative).

[tout, qout] = RK2(qcFunc, t0, qc0, h, a, tf); %Numerical method 

qcExact = ( exp(-tout/(R*C)) .* ( C^2 * w^2 * R^2 + 2000000 * C^2 * w * R * V * exp(tout/(R*C)) .* sin(w*tout) + 2000000 * C * V * exp(tout/(R*C)) .* cos(w*tout) - 2000000 * C * V + 1 ) ) ./ (2000000 * (C^2 * w^2 * R^2 + 1));        
% MY GOD this exact solution was hard. I better get the nobel prize for
% this.

error = abs(qout - qcExact); % get the error between the numerical and the exact one.

plot(tout, qout, 'b'); %plot the numerical graph and make sure it comes up and 'holds' so we can compare.
ylabel("qc(t)");
xlabel("t");
title("numerical");
figure;

plot(tout, qcExact, 'r'); %now plot the exact graph
ylabel("qc(t)");
xlabel("t");
title("exact");
figure;

plot(tout, error); %now plot the error 
ylabel("error");
xlabel("t");
title("error");
figure;

loglog(tout, error); %now plot the error using a log-log plot to see that mag.
ylabel("log(error)");
xlabel("log(t)");

title("loglog");









