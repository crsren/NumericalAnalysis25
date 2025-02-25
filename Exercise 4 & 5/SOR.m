% Successive over-relaxation (used by testbench.h)

function [iterations] = SOR(a,h,n,m,Ui_min,Ui_max,Umin_j,Umax_j,g)
% Inputs:
% SOR( Relaxation factor α, grid square size h, 2 grid dimensions
% 4 boundary conditions, 0 )
    
% Initialize grid (all zero)
x = 0:h:n;
y = 0:h:m;
U = zeros(n+1,m+1);

% Set boundary conditions
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
% Define timeout (maximum iterations)
timeout = 10000;

iterations = -1; % counter
precise = false;

% Apply relexation method on all interior points until the residue r < Ɛ
while ~precise
    precise = true;
    
    for i = 2 : length(x)-1
        for j = 2 : length(y)-1
            
            % Calculate residue r
            r = ( U(i+1,j) + U(i-1,j) + U(i,j+1) + U(i,j-1) - 4*U(i,j)  - (h^2)*g(i,j) ) / (4);
            U(i,j) = U(i,j) + a*r;
            
            if abs(r) >= epsilon
                precise = false;
            end
            
        end
    end
    
    iterations = iterations+1;
    if(iterations >= timeout)
        disp("Time out.");
        break;
    end
end

disp("Done after " + iterations + " iterations (α = " + a + ")");

% Plotting
[X,Y] = meshgrid(x,y);
figure; surf(X,Y,U); % 3D
% figure; contour(X,Y,U); % from above
% figure; meshc(X,Y,U) %mesh: contour + surf
xlabel('x'); ylabel('y'); zlabel('U(x,y)'); title("α = " + a);
