%Setting up the T domain
h = 0.1;
tf = 2;
t0 = 0;
N = tf/h;
%Setting up Variables
qc = zeros(1,N);
qc(1) = 500*10^(-9);
qc_grad = zeros(1,N);
Vout = zeros(1,N);
t = zeros(1,N);
t(1) = t0;
F = @f2;

[t, qc, qc_grad]= RK4(F, h, t, qc, qc_grad);

Vout = R*qc_grad;

plot(t, Vout);
xlabel('Time');
ylabel('Potential Difference');

function [yout] = f2(time, qc, qc_grad)
    R = 250;
	C = 3.5*10^(-6);
    L = 600*10^(-3);
    yout = (1/L)*(Vin(time) - R*qc_grad - (1/C)*qc);
end

function [volt_out] = Vin(time);
    volt_out = 5;
end


