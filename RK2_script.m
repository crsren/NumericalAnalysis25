y0 = 5;
x0 = 0;
h = 0.00001;
xf = 0.01;
a = 0.5;
%a = 1; 
vi = @(x) 2.5 * stepfun(x, 0);
%vi = @(x) 2.5 * exp(-x^2/0.0001);
%vi = @(x) 2.5 * exp(-x/0.0001);
%T =0.00001;vi = @(x) 2.5 * sin(2*pi*x/T);
f = @(x, y) 10000 * (vi(x) - y);
[y, x] = RK2(x0, y0, h, xf, a, f);
plot(x, y)