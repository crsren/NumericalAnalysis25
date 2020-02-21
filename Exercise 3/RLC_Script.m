R = 250;
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
F1 = @f1;
F2 = @f2;
%Input Signal
Vin = 5;


[qc, qc_grad]= RK4(F1, F2, h, t, qc, qc_grad);

for i=1: (length(t)-1)
    t(i+1) = t(i)+h;
    Vout(i) = R*abs(qc_grad(i));
end

plot(t, Vout);
set(gca, 'YScale', 'log')
xlabel('Time');
ylabel('Amplitude');

   
function [yout] = f1(Vin, qc, qc_grad)
    yout = qc_grad;
end

function [yout] = f2(Vin, qc, qc_grad)
    R = 250;
	C = 3.5*10^(-6);
    L = 600*10^(-3);
    yout = (1/L)*(Vin - R*qc_grad - (1/C)*qc);
end

