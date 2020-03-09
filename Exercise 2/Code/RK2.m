function[x, y] = RK2(func, x0, y0, h, a, xf)
    
    a = a; % determining the constants the determines the RK method. all based on a.
    b = 1 - a; % determining the constants the determines the RK method. all based on a.
    p = 0.5/b; % determining the constants the determines the RK method. all based on a.
    q = 0.5/b; % determining the constants the determines the RK method. all based on a.
    
    N = ceil((xf-x0)/h); %number of steps
    
    x = zeros(1, N); %create array populated with zeros. 1xN
    y = zeros(1, N); %create array populated with zeros. 1xN
    
    x(1) = x0; %Initialize first value to t=0
    y(1) = y0; %Initialize first value to t=0
    
    for i = 1:N-1 %iteration fun
        
        k1 = func(x(i), y(i)); %determine k1
        k2 = func( x(i)+h, y(i) + k1*h ); %determine k2
        
        x(i+1) = x(i) + h; %next x is obviously previous x + the step distance (h).
        y(i+1) = y(i) + h * (a*k1 + b*k2); %RK definition.
            
    end %end loop
end %end function