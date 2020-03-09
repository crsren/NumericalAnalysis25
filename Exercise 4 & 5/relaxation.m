% Normal relaxation method
% later: function [iterations,U] = relaxation(h,n,m,epsilon,Ui_min,Ui_max,Umin_j,Umax_j)
  
% [ ] examples for other poisson functions ( = g(x) instead of 0)
% [ ] understand SOR

close all ; clc

% d2u/dx2 + d2u/dy2 = g(x,y) → Poisson Equation
% Case g(x,y) = 0 → Laplace Equation
g = @(x,y) 0; 

% Initialize grid (all zero)
h = 0.1; n = 20; m = 20;
x = 0:h:h*n;
y = 0:h:h*m;
U = zeros(n+1,m+1);

% Set boundary conditions
Ui_min = @(x) x;
Ui_max = @(x) x;
Umin_j = @(y) y;
Umax_j = @(y) y;

for i = 1 : length(x)
    U(i ,1) = Ui_min(i);
    U(i,length(y)) = Ui_max(i);
end
for j = 1 : length(y)
    U(1, j) = Umin_j(j);
    U(length(x), j) = Umax_j(j);
end

% Define termination condition Ɛ (maximum error)
epsilon = 1e-20;

% Apply relexation method on all interior points until the residue r < Ɛ
iterations = -1; % counter
precise = false;

while ~precise
    precise = true;
    
    for i = 2 : length(x)-1
        for j = 2 : length(y)-1
            
            % Calculate residue r
            r = ( U(i+1,j) + U(i-1,j) + U(i,j+1) + U(i,j-1) - 4*U(i,j) - (h^2)*g(i,j) ) / (4);
            U(i,j) = U(i,j) + r;
            
            if abs(r) >= epsilon
                %disp(r);
                precise = false;
            end
            
        end
    end
    
    iterations = iterations+1;
end

disp("Done after " + iterations + " iterations.");


% Plotting
[X,Y] = meshgrid(x,y);
figure; surf(X,Y,U); % 3D
xlabel('x'); ylabel('y'); zlabel('Potential in Volts');

figure; contour(X,Y,U); % from above
%figure; meshc(X,Y,U) %mesh: contour + surf
