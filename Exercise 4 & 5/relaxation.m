close all ; clc

% d2u/dx2 + d2u/dy2 = g(x,y) → Poisson Equation
g = @(x,y) 0;

% Initialize grid (all zero)
h = 0.01;n = 1; m = 1;
x = 0:h:n;
y = 0:h:m;
U = zeros(length(x),length(y));

% Set boundary conditions
Ui_min = @(x) x;
Ui_max = @(x) x;
Umin_j = @(y) y;
Umax_j = @(y) y;

for i = 1 : length(x)
    U(i ,1) = Ui_min(x(i));
    U(i,length(y)) = Ui_max(x(i));
end
for j = 1 : length(y)
    U(1, j) = Umin_j(y(j));
    U(length(x), j) = Umax_j(y(j));
end

% Define termination condition Ɛ (maximum error)
epsilon = 1e-4;

iterations = -1; % counter
precise = false;

% Apply relexation method on all interior points until the residue r < Ɛ
while ~precise
    precise = true;
    
    for i = 2 : length(x)-1
        for j = 2 : length(y)-1
            
            % Calculate residue r
            r = ( U(i+1,j) + U(i-1,j) + U(i,j+1) + U(i,j-1) - 4*U(i,j)  - (h^2)*g(i,j) ) / (4);
            U(i,j) = U(i,j) + r;
            
            if abs(r) >= epsilon
                precise = false;
            end
            
        end
    end
    
    iterations = iterations+1;
end

disp("Done after " + iterations + " iterations with " + (length(x)-2)^2 + " calculations each.");

% Plotting
[X,Y] = meshgrid(x,y);
figure; surf(X,Y,U); % 3D
%figure; meshc(X,Y,U) %mesh: contour + surf
xlabel('x'); ylabel('y'); zlabel('U(x,y)');
title("Resolution h = 0.01");
view(-70,20); saveas(gcf,'h001.png');