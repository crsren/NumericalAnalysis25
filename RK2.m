function [y, x] = RK2(x0, y0, h, xf, a, f)
   
    c1 = a;
    c2 = 1 - a;
    u = 1/(2 * a);
    m = 1/(2 * a);
    N = ceil((xf - x0)/h);
    y = zeros(1, N);
    y(1) = y0;
    x = zeros(1, N);
    x(1) = x0;
    for i = 1:N - 1
        K1 = f(x(i), y(i));
        K2 = f(x(i) + m * h, y(i) + u * h * K1);
        y(i + 1) = y(i) + h * (c1 * K1 + c2 * K2);
        x(i + 1) = x(i) + h;
    end
end