% Normal relaxation method
% later: function [iterations,U] = relaxation(h,n,m,epsilon,Ui_min,Ui_max,Umin_j,Umax_j)

%close all ; clc

% d2u/dx2 + d2u/dy2 = g(x,y) → Poisson Equation
% Case g(x,y) = 0 → Laplace Equation
g = @(x,y) 0; % !! Poisson eqn other than Laplace not working yet !!

% Initialize grid (all zero)
h = 0.02; n = 1; m = 1;
x = 0:h:n;
y = 0:h:m;
U = zeros(n+1,m+1);


% Set boundary conditions
Ui_min = @(x) discon(x);
Ui_max = @(x) discon(x);
Umin_j = @(y) discon(y);
Umax_j = @(y) discon(y);

for i = 1 : length(x)
    U(i ,1) = Ui_min(x(i));
    U(i,length(y)) = Ui_max(x(i));
end
for j = 1 : length(y)
    U(1, j) = Umin_j(y(j));
    U(length(x), j) = Umax_j(y(j));
end

iterations = (1:1:5000);
% r1map = zeros(length(iterations),1);
% r2map = zeros(length(iterations),1);
%r3map = zeros(length(iterations),1);
rdist = zeros(n+1,m+1);

for n=1:length(iterations)
    
    for i = 2 : length(x)-1
        for j = 2 : length(y)-1
            
            % Calculate residue r
            r = ( U(i+1,j) + U(i-1,j) + U(i,j+1) + U(i,j-1) - 4*U(i,j)  - (h^2)*g(i,j) ) / (4);
            U(i,j) = U(i,j) + r;
            
            if n==length(iterations)
                rdist(i,j) = r;
            end
            
%             if ((i == 99) && (j == 99))
%                 r1map(n,1) = r;
%             end
%             if ((i == 76) && (j == 76))
%                 r2map(n,1) = r;
%             end
%             if ((i == 51) && (j == 51))
%                 r3map(n,1) = r;
%             end
            
        end
    end
end

disp("Done after " + iterations + " iterations.");


% Plotting
[X,Y] = meshgrid(x,y);
figure; surf(X,Y,U); % 3D
xlabel('x'); ylabel('y'); zlabel('U(x,y)');
%saveas(gcf,'RDgraphsin5000it.png');

[A,B] = meshgrid(h:h:1,h:h:1);
figure; surf(A,B,rdist);
xlabel('x'); ylabel('y'); zlabel('Residue r');
title("Residue distribution after " + length(iterations) + " iterations");
%saveas(gcf,'RDsin5000it.png');

% figure; hold on;
% plot(iterations, r1map);
% plot(iterations, r2map);
% plot(iterations, r3map);
% e = ones(size(iterations))*0.015;
% plot(iterations, e, "--black");
% ylim([0 0.15]);
% title("")
% xlabel("Iterations");
% ylabel("Residue r");
%legend("i=j= 5","i=j= 75","i=j=100","Ɛ");

%figure; contour(X,Y,U); % from above
%figure; meshc(X,Y,U) %mesh: contour + surf

%Fn for special boundary condition with non-corner discontinuity
function [yy] = discon(xx)
    if xx < 0.25
        yy = 1;
    elseif xx < 0.5
        yy = 0;
    elseif xx < 0.75
        yy = 1;
    else
        yy = 0;
    end
end
