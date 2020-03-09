% Successive over-relaxation

function [iterations] = SOR(a,h,n,m,Ui_min,Ui_max,Umin_j,Umax_j,g) % input: alpha weighting of residue
    
% Initialize grid (all zero)
x = 0:h:h*n;
y = 0:h:h*m;
U = zeros(n+1,m+1);


for i = 1 : length(x)
    U(i ,1) = Ui_min(x(i));
    U(i,length(y)) = Ui_max(x(i));
end
for j = 1 : length(y)
    U(1, j) = Umin_j(y(j));
    U(length(x), j) = Umax_j(y(j));
end

% Define termination condition Ɛ (maximum error)
epsilon = 1e-2;
% Define timeout (maximum iterations)
timeout = 10000;

% Apply relexation method on all interior points until the residue r < Ɛ
iterations = -1; % counter
precise = false;

while ~precise
    precise = true;
    
    for i = 2 : length(x)-1
        for j = 2 : length(y)-1
            
            % Calculate residue r
            r = ( U(i+1,j) + U(i-1,j) + U(i,j+1) + U(i,j-1) - 4*U(i,j)  - (h^2)*g(i,j) ) / (4);
            U(i,j) = U(i,j) + a*r;
            
            if abs(r) >= epsilon
                %disp(r);
                precise = false;
            end
            
        end
    end
    
    iterations = iterations+1;
    if(iterations >= timeout)
        disp("Time out.");
        break;
    end
    %disp(iterations)
end

disp("Done after " + iterations + " iterations.");


% % Plotting
% [X,Y] = meshgrid(x,y);
figure; surf(X,Y,U); % 3D
% xlabel('x'); ylabel('y'); zlabel('Potential in Volts');
% 
% figure; contour(X,Y,U); % from above
% %figure; meshc(X,Y,U) %mesh: contour + surf
