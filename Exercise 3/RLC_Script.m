R = 250;
L = 600*10^(-3);
C = 3.5*10^(-6);

qc0 = 500*10^(-9);

%System properties
h = 0.00001;
tf = 0.05;
t0 = 0;
N = tf/h;

%Setting up Variables
qc = zeros(1,N);
qc(1) = qc0;
qc_grad = zeros(1,N);
Vout = zeros(1,N);
t = zeros(1,N);
t(1) = t0;

F = @y;

%Call to RK4
[t, qc, qc_grad] = RK4(F, h, t, qc, qc_grad);

Vout = R*qc_grad;

%obtaining plot
plot(t, Vout);
xlabel('Time (s)');
ylabel('Potential Differencev (V)');
hold on;




%System function. other one doesn't have to be defined (just qc_dash)
%change number after Vin for different imput signals shown below
function [yout] = y(t, qc, qc_grad)
    R = 250;
    L = 600*10^(-3);
    C = 3.5*10^(-6);
    yout = (Vin4(t) - (qc/C) - (R*qc_grad))/L;
end

% (1) heaviside Vin
function [volt_out] = Vin1(time)
    volt_out = 5;
end


% (2) impulsive signal
function [volt_out] = Vin2(time)
    volt_out = 5*exp(-(time^2)/3);
end

% (3) square wave 
function [volt_out] = Vin3(time)
    freq = 109;
    volt_out = 5*square(2*pi*freq*time);
end

% (4) sine wave
function [volt_out] = Vin4(time)
    freq = 109;
    volt_out = 5*sin(2*pi*freq*time);
end


