% Normal relaxation method
% later: function [iterations,U] = relaxation(h,n,m,epsilon,Ui_min,Ui_max,Umin_j,Umax_j)
  
%-> Next: implement functions as boundary conditions

close all ; clc

% Initialize grid (all zero)
h = 1; n = 30; m = 30;
x = 0:h:n;
y = 0:h:m;
U = zeros(n+1,m+1);

% Set boundary conditions
for i = 1 : length(x)
    U(i ,1) = 20;
    U(i,length(y)) = 20;
end
for j = 1 : length(y)
    U(1, j) = 10;
    U(length(x), j) = 0;
end

% Define termination condition Ɛ (maximum error)
epsilon = 0.01;

% Apply relexation method on all interior points until the residue r < Ɛ
iterations = -1; % counter
precise = false;

while ~precise
    precise = true;
    
    for i = 2 : length(x)-1
        for j = 2 : length(y)-1
            
            % Calculate residue r
            r = ( U(i+1,j) + U(i-1,j) + U(i,j+1) + U(i,j-1) - 4*U(i,j) ) / (4*h);
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
figure; contour(X,Y,U); % from above

%subplot(2,2,3); contour3(X,Y,U) %shows a contour plot with colored lines xlabel('x'); ylabel('y'); zlabel('Potential in Volts');
%subplot(2,2,4); meshc(X,Y,V) %shows a mesh plot of the solution xlabel('x'); ylabel('y'); zlabel('Potential in Volts');
