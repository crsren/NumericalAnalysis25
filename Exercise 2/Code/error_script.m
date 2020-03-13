clear;
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
      a = 0.3;
%*********************************************



Vi = @(t) V * cos(w*t); %input cosine wave with given values.

qcFunc = @(t, qc) ( Vi(t)/R ) - ( 1/(R*C) ) * qc; %we get the d/dt of qc (derivative).

[tout, qout] = RK2(qcFunc, t0, qc0, h, a, tf); %Numerical method 

qcExact = ( exp(-tout/(R*C)) .* ( C^2 * w^2 * R^2 + 2000000 * C^2 * w * R * V * exp(tout/(R*C)) .* sin(w*tout) + 2000000 * C * V * exp(tout/(R*C)) .* cos(w*tout) - 2000000 * C * V + 1 ) ) ./ (2000000 * (C^2 * w^2 * R^2 + 1));        
% MY GOD this exact solution was hard. I better get the nobel prize for
% this.

error = qcExact - qout; % get the error between the numerical and the exact one.

plot(tout, qout, 'b'); %plot the numerical graph and make sure it comes up and 'holds' so we can compare.
ylabel("qc(t)");
xlabel("t");
title("numerical solution");

figure;
plot(tout, qcExact, 'r'); %now plot the exact graph
ylabel("qc(t)");
xlabel("t");
title("exact solution");
figure;

plot(tout, error); %now plot the error 
ylabel("error");
xlabel("t");
title("error");
figure;

%%%%     For this part. These first three for-loops are for doing the
%%%%     log-log plots of the errors vs step-size. The last for-loop is for
%%%%     doing a semi-log where the error is not logarithmized but the
%%%%     step-size is for convenience. There are some random 'a = xx'.
%%%%     Disregard these as only I know the right combination.
%%%%     Heun's method Error Loglog

%  a = 0.5;
% 
% for i = 16:25
%     
%     step = 2^(-i); %Trying different step sizes (h).
%     [etout, eqout] = RK2(qcFunc, t0, qc0, step, a, tf); %RK2 functions with different step size h. added extra e for fun and syntax.
%     qcExact = ( exp(-etout/(R*C)) .* ( C^2 * w^2 * R^2 + 2000000 * C^2 * w * R * V * exp(etout/(R*C)) .* sin(w*etout) + 2000000 * C * V * exp(etout/(R*C)) .* cos(w*etout) - 2000000 * C * V + 1 ) ) ./ (2000000 * (C^2 * w^2 * R^2 + 1));
%     error = max(abs(qcExact - eqout)); %redefinition doesn't seem to bother matlab. Also had to call exact function again.
%     plot(log(step), log(error), 'b*'); % doesn't seem to make a big difference doing log-log vs plot log.
%     ylabel("log Error = max|Exact - Numerical|"); % add some labels.
%     xlabel("log step-size");
%     %Uncomment the title for when you make the compare graph
%     title("Comparing Heun's method log(error) (red) with Midpoint method's (blue) with Random method's (green)");
%     %title("Heun's method. Log-Log graph of error to step-size");
%     hold on; %10 iterations of plotting
%     
% end

%%%%     Midpoint method Error Loglog

% a = 0.0;
% 
% for i = 16:25
% 
%     step = 2^(-i); %Trying different step sizes (h).
%     [etout, eqout] = RK2(qcFunc, t0, qc0, step, a, tf); %RK2 functions with different step size h. added extra e for fun and syntax.
%     qcExact = ( exp(-etout/(R*C)) .* ( C^2 * w^2 * R^2 + 2000000 * C^2 * w * R * V * exp(etout/(R*C)) .* sin(w*etout) + 2000000 * C * V * exp(etout/(R*C)) .* cos(w*etout) - 2000000 * C * V + 1 ) ) ./ (2000000 * (C^2 * w^2 * R^2 + 1));
%     error = max(abs(qcExact - eqout)); %redefinition doesn't seem to bother matlab. Also had to call exact function again.
%     plot(log(step), log(error), 'r*'); % doesn't seem to make a big difference doing log-log vs plot log.
%     ylabel("log Error = max|Exact - Numerical|"); % add some labels.
%     xlabel("log step-size");
%     %Uncomment the title for when you make the compare graph if needed.
%     title("Comparing Heun's method log(error) (red) with Midpoint method's (blue) with Random method's (green)");
%   % title("Midpoint method. Log-Log graph of error to step-size");
%     hold on; %10 iterations of plotting
%     
% end


%%%%     Random method Error Loglog

% a = 0.3;
% 
% for i = 16:25
% 
%     step = 2^(-i); %Trying different step sizes (h).
%     [etout, eqout] = RK2(qcFunc, t0, qc0, step, a, tf); %RK2 functions with different step size h. added extra e for fun and syntax.
%     qcExact = ( exp(-etout/(R*C)) .* ( C^2 * w^2 * R^2 + 2000000 * C^2 * w * R * V * exp(etout/(R*C)) .* sin(w*etout) + 2000000 * C * V * exp(etout/(R*C)) .* cos(w*etout) - 2000000 * C * V + 1 ) ) ./ (2000000 * (C^2 * w^2 * R^2 + 1));
%     error = max(abs(qcExact - eqout)); %redefinition doesn't seem to bother matlab. Also had to call exact function again.
%     plot(log(step), log(error), 'g*'); % doesn't seem to make a big difference doing log-log vs plot log.
%     ylabel("log Error = max|Exact - Numerical|"); % add some labels. 
%     xlabel("log step-size h");
%     %Uncomment the title for when you make the compare graph
%     title("Comparing Heun's method log(error) (red) with Midpoint method's (blue) with Random method's (green)");
%     %title("Random method. Log-Log graph of error to step-size");
%     hold on; %10 iterations of plotting
% 
% end

%%% All three methods for the error vs step-size (NOT LOGLOG! but semilog)



for i = 16:25

    step = 2^(-i); %Trying different step sizes (h).
    [etout, eqout] = RK2(qcFunc, t0, qc0, step, a, tf); %RK2 functions with different step size h. added extra e for fun and syntax.
    qcExact = ( exp(-etout/(R*C)) .* ( C^2 * w^2 * R^2 + 2000000 * C^2 * w * R * V * exp(etout/(R*C)) .* sin(w*etout) + 2000000 * C * V * exp(etout/(R*C)) .* cos(w*etout) - 2000000 * C * V + 1 ) ) ./ (2000000 * (C^2 * w^2 * R^2 + 1));
    error = max(abs(qcExact - eqout)); %redefinition doesn't seem to bother matlab. Also had to call exact function again.
     plot(log(step), error, 'b*');
     ylabel("Error = max|Exact - Numerical|");
     xlabel("step-size h");
%     title("Heun's method error vs log(step-size)");
%     title("Midpoint method error vs log(step-size)");
     title("Random method error vs log(step-size)");
     hold on;
end










