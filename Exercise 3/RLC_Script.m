%System properties
h = 2^(-15);
tf = 0.05;
t0 = 0;
N = tf/h;
qc0 = 500*10^(-9);

%Setting up Variables
qc(1) = qc0;
qc_grad(1) = 0;
Vout(1) = 0;
t(1) = t0;

F1 = @y;
F2 = @z;

%Call to RK4
[Vout, t, qc, qc_grad] = RK4(Vout, F1, F2, h, t, qc, qc_grad, N);

%obtaining plot
plot(t, Vout);
xlabel('Time (s)');
ylabel('Potential Differencev (V)');
title('Step signal with Vin = 5V')
hold on;


%System function.
function [yout] = z(t, qc, qc_grad)
    yout = qc_grad;
end

%change number after Vin for different imput signals shown below
function [yout] = y(t, qc, qc_grad)
    R = 250;
    L = 600*10^(-3);
    C = 3.5*10^(-6);
    yout = (Vin1(t) - (qc/C) - (R*qc_grad))/L;
end

% (1) heaviside Vin
function [volt_out] = Vin1(time)
    volt_out = 5*heaviside(time);
end


% (2) impulsive signal
function [volt_out] = Vin2(time)
    volt_out = 5*exp(-(time^2)/(3*10^(-6)));
end

% (3) square wave 
function [volt_out] = Vin3(time)
    freq = 500;
    volt_out = 5*square(2*pi*freq*time);
end

% (4) sine wave
function [volt_out] = Vin4(time)
    freq = 5;
    volt_out = 5*sin(2*pi*freq*time);
end


