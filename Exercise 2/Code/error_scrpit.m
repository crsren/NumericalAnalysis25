clear;
R = 1e3; %resistance.
C = 100 * 1e-9; %capacitance.
qc0 = 500 * 1e-9; %charge at t0.
t0 = 0; %value of time when time is zero (obviously).
h = 0.00001; %step size.
tf = 0.001; %interval [0, 0.01].


%****************Heun's Method****************
a = 0.5;
%*********************************************

%***************Midpoint Method***************
% a = 0.0;
%*********************************************

%*****************Random RK2******************
% a = 0.3;
%*********************************************


%*****************Sinusoidal Input Signal******************
V = 5.0; %amplitude (voltage).
T = 100 * 1e-6; %period of T (1/T = frequency). 2pi*frequency = rotationshstighet (angular frequency?).
freq = 1/T; % frequency 
w = 2 * pi * freq; %angular frequency
Vi = @(t) V * cos(w*t); %input cosine wave with given values.
%**********************************************************

%*****************Decay Input Signal******************
% V = 2.5;
% tau = 100 * 1e-6;
% Vi = @(t) V * exp(-t/tau);
%*****************************************************

qcFunc = @(t, qc) ( Vi(t)/R ) - ( 1/(R*C) ) * qc; %we get the d/dt of qc (derivative).

[tout, qout] = RK2(qcFunc, t0, qc0, h, a, tf); %Numerical method

%*****************Exact Solutions of the 2 Input Signal******************
qcExact = ( exp(-tout/(R*C)) .* ( C^2 * w^2 * R^2 + 2000000 * C^2 * w * R * V * exp(tout/(R*C)) .* sin(w*tout) + 2000000 * C * V * exp(tout/(R*C)) .* cos(w*tout) - 2000000 * C * V + 1 ) ) ./ (2000000 * (C^2 * w^2 * R^2 + 1));
% qcExact = V * tout .* exp(-tout/(R*C)) / R + 2 * V * C * exp(-tout/(R*C));
%************************************************************************

error = qcExact - qout; % get the error between the numerical and the exact one.

plot(tout, qout, 'b'); %plot the numerical graph and make sure it comes up and 'holds' so we can compare.
hold on;
plot(tout, qcExact, 'r'); %now plot the exact graph
hold on;
plot(tout, error, 'g'); %now plot the error 
ylabel("qc(t)");
xlabel("t");
title("Heun's Method: Error Analysis");
% title("Midpoint Method: Error Analysis");
% title("Random Method: Error Analysis");
legend("Numerical Solution", "Exact Solution", "Approximation Error");
figure;

plot(tout, error, 'g'); %now plot the error 
ylabel("Error");
xlabel("t");
title("Heun's Method: Error Analysis");
% title("Midpoint Method: Error Analysis");
% title("Random Method: Error Analysis");
figure;



%%%%     For this part. These first three for-loops are for doing the
%%%%     log-log plots of the errors vs step-size. The last for-loop is for
%%%%     doing a semi-log where the error is not logarithmized but the
%%%%     step-size is for convenience. There are some random 'a = xx'.
%%%%     Disregard these as only I know the right combination.


%%%%     Heun's method Error Loglog
 a = 0.5;

for i = 16:25
    
    step = 2^(-i); %Trying different step sizes (h).
    [etout, eqout] = RK2(qcFunc, t0, qc0, step, a, tf);
    
    %*****************Exact Solutions of the 2 Input Signal******************
    qcExact = ( exp(-etout/(R*C)) .* ( C^2 * w^2 * R^2 + 2000000 * C^2 * w * R * V * exp(etout/(R*C)) .* sin(w*etout) + 2000000 * C * V * exp(etout/(R*C)) .* cos(w*etout) - 2000000 * C * V + 1 ) ) ./ (2000000 * (C^2 * w^2 * R^2 + 1));
%     qcExact = V * etout .* exp(-etout/(R*C)) / R + 2 * V * C * exp(-etout/(R*C));
    %************************************************************************

    error = max(abs(qcExact - eqout));
    p1 = plot(log(step), log(error), 'b*'); % doesn't seem to make a big difference doing log-log vs plot log.
    hold on; %10 iterations of plotting
    
end
hold on;

%%%     Midpoint method Error Loglog
a = 0.0;

for i = 16:25

    step = 2^(-i); %Trying different step sizes (h).
    [etout, eqout] = RK2(qcFunc, t0, qc0, step, a, tf);
    
    %*****************Exact Solutions of the 2 Input Signal******************
    qcExact = ( exp(-etout/(R*C)) .* ( C^2 * w^2 * R^2 + 2000000 * C^2 * w * R * V * exp(etout/(R*C)) .* sin(w*etout) + 2000000 * C * V * exp(etout/(R*C)) .* cos(w*etout) - 2000000 * C * V + 1 ) ) ./ (2000000 * (C^2 * w^2 * R^2 + 1));
%     qcExact = V * etout .* exp(-etout/(R*C)) / R + 2 * V * C * exp(-etout/(R*C));
    %************************************************************************

    error = max(abs(qcExact - eqout));
    p2 = plot(log(step), log(error), 'r*'); % doesn't seem to make a big difference doing log-log vs plot log.
    hold on; %10 iterations of plotting
    
end
hold on;

%%%     Random method Error Loglog
a = 0.3;

for i = 16:25

    step = 2^(-i); %Trying different step sizes (h).
    [etout, eqout] = RK2(qcFunc, t0, qc0, step, a, tf);
    
    %*****************Exact Solutions of the 2 Input Signal******************
    qcExact = ( exp(-etout/(R*C)) .* ( C^2 * w^2 * R^2 + 2000000 * C^2 * w * R * V * exp(etout/(R*C)) .* sin(w*etout) + 2000000 * C * V * exp(etout/(R*C)) .* cos(w*etout) - 2000000 * C * V + 1 ) ) ./ (2000000 * (C^2 * w^2 * R^2 + 1));
%     qcExact = V * etout .* exp(-etout/(R*C)) / R + 2 * V * C * exp(-etout/(R*C));
    %************************************************************************

    error = max(abs(qcExact - eqout));
    p3 = plot(log(step), log(error), 'g*'); % doesn't seem to make a big difference doing log-log vs plot log.
    hold on; %10 iterations of plotting

end

ylabel("Log (Error)"); % add some labels.
xlabel("log (Step-Size)");
legend([p1, p2, p3], "Heun's method", "Midpoint method", "Random method", 'Location','northwest');
title("Log (Error) against Log (Step-Size)");
figure;



% %% All three methods for the error vs step-size (NOT LOGLOG! but semilog)
%%%%     Heun's method Error SemiLog
a = 0.5;
for i = 16:25

    step = 2^(-i); %Trying different step sizes (h).
    [etout, eqout] = RK2(qcFunc, t0, qc0, step, a, tf);
    
    %*****************Exact Solutions of the 2 Input Signal******************
    qcExact = ( exp(-etout/(R*C)) .* ( C^2 * w^2 * R^2 + 2000000 * C^2 * w * R * V * exp(etout/(R*C)) .* sin(w*etout) + 2000000 * C * V * exp(etout/(R*C)) .* cos(w*etout) - 2000000 * C * V + 1 ) ) ./ (2000000 * (C^2 * w^2 * R^2 + 1));
%     qcExact = V * etout .* exp(-etout/(R*C)) / R + 2 * V * C * exp(-etout/(R*C));
    %************************************************************************

    error = max(abs(qcExact - eqout));
    p4 = plot(log(step), error, 'b*');
    hold on;
end
hold on;

%%%%     Midpoint method Error SemiLog
a = 0.0;
for i = 16:25

    step = 2^(-i); %Trying different step sizes (h).
    [etout, eqout] = RK2(qcFunc, t0, qc0, step, a, tf);
    
    %*****************Exact Solutions of the 2 Input Signal******************
    qcExact = ( exp(-etout/(R*C)) .* ( C^2 * w^2 * R^2 + 2000000 * C^2 * w * R * V * exp(etout/(R*C)) .* sin(w*etout) + 2000000 * C * V * exp(etout/(R*C)) .* cos(w*etout) - 2000000 * C * V + 1 ) ) ./ (2000000 * (C^2 * w^2 * R^2 + 1));
%     qcExact = V * etout .* exp(-etout/(R*C)) / R + 2 * V * C * exp(-etout/(R*C));
    %************************************************************************

    error = max(abs(qcExact - eqout));
    p5 = plot(log(step), error, 'r*');
    hold on;
end
hold on;

%%%%     Random method Error SemiLog
a = 0.3;
for i = 16:25

    step = 2^(-i); %Trying different step sizes (h).
    [etout, eqout] = RK2(qcFunc, t0, qc0, step, a, tf);
    
    %*****************Exact Solutions of the 2 Input Signal******************
    qcExact = ( exp(-etout/(R*C)) .* ( C^2 * w^2 * R^2 + 2000000 * C^2 * w * R * V * exp(etout/(R*C)) .* sin(w*etout) + 2000000 * C * V * exp(etout/(R*C)) .* cos(w*etout) - 2000000 * C * V + 1 ) ) ./ (2000000 * (C^2 * w^2 * R^2 + 1));
%     qcExact = V * etout .* exp(-etout/(R*C)) / R + 2 * V * C * exp(-etout/(R*C));
    %************************************************************************

    error = max(abs(qcExact - eqout));
    p6 = plot(log(step), error, 'g*');
    hold on;
end

ylabel("Error");
xlabel("Step-Size");
legend([p4, p5, p6], "Heun's method", "Midpoint method", "Random method", 'Location','northwest');
title("Error against Log (Step-Size)");