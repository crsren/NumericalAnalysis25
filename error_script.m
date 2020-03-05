y0 = 5;
x0 = 0;
h = 0.000001;
xf = 0.01;
a = 1; 
T = 0.0001;
vi = @(x) 5 * cos(2*pi*x/T);
f = @(x, y) 10000 * (vi(x) - y);
[y, x] = RK2(x0, y0, h, xf, a, f);
y_exact = (5 * cos(2*pi*x/T) + 10 * pi * sin(2*pi*x/T) + 20 * pi * pi * exp(-10000.* x))/(1 + 4 * pi *pi);
error = y - y_exact;
figure
plot(x, y);
title('solution');
figure
plot(x, y_exact);
title('exact');
figure
plot(x, error);
title('error');