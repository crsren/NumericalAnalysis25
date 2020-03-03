function [xout, yout] = RK2(func, x_i, y_i, x_f, h)
    
    N = round( (x_f - x_i)/h ); % use this determine size of arrays

    %****************Heun's Method****************
    a = 0.5;
    %*********************************************

    %***************Midpoint Method***************
    % a = 0;
    %*********************************************

    %*****************Random RK2*****************
    % a = 0.3;
    %*********************************************

    b = 1 - a;
    p = 0.5/b;
    q = p;

    yout = zeros(1,N);
    xout = zeros(1,N); %set up arrays
    xout(1) = x_i;
    yout(1) = y_i; % initialize arrays

    for j = 1:N-1 % loop for N-1 steps
        k1 = feval( func, xout(j), yout(j) );
        k2 = feval( func, (xout(j) + p*h), (yout(j) + q*k1*h) );
        yout(j+1) = yout(j) + h*(a*k1 + b*k2); % next value of y calculated from previous values of x,y
        xout(j+1) = xout(j) + h; % increase x by stepsize
    end
end