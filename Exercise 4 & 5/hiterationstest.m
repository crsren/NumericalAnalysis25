% Normal relaxation method
% later: function [iterations,U] = relaxation(h,n,m,epsilon,Ui_min,Ui_max,Umin_j,Umax_j)

%close all ; clc

% d2u/dx2 + d2u/dy2 = g(x,y) → Poisson Equation
% Case g(x,y) = 0 → Laplace Equation
g = @(x,y) 0; % !! Poisson eqn other than Laplace not working yet !!

% Initialize grid (all zero)

% Set boundary conditions
Ui_min = @(x) x;%(x-100)^2;
Ui_max = @(x) x;%(x-100)^2;
Umin_j = @(y) y;%(y-100)^2;
Umax_j = @(y) y;%(y-100)^2;

hset = 0.005:0.005:0.175;
recordings = zeros(1, length(hset));
recordingsi = recordings;
recordingsc = recordings;

for h=1:length(hset)
    
    n = 1; m = 1;
    x = 0:hset(h):n;
    y = 0:hset(h):m;
    U = zeros(length(x),length(y));

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
    
    % Apply relexation method on all interior points until the residue r < Ɛ
    iterations = -1; % counter
    precise = false;
    
    while ~precise
        precise = true;
        
        for i = 2 : length(x)-1
            for j = 2 : length(y)-1
                
                % Calculate residue r
                r = ( U(i+1,j) + U(i-1,j) + U(i,j+1) + U(i,j-1) - 4*U(i,j)  - (h^2)*g(i,j) ) / (4);
                U(i,j) = U(i,j) + r;
                
                if abs(r) >= epsilon
                    %disp(r);
                    precise = false;
                end
                
            end
        end
        
        iterations = iterations+1;
        %disp(iterations)
    end
    
    disp(hset(h) + ": Done after " + iterations + " iterations with " + (length(x)-2)^2 + " calculations each.");
    
    recordings(h) = iterations*(length(x)-2)^2;
    recordingsi(h) = iterations;
    recordingsc(h) = (length(x)-2)^2;
end
figure; plot(hset, recordings, "b");
xlabel("Resolution h"); ylabel("Total calculations");
figure; plot(hset, recordingsi, "r");
xlabel("Resolution h"); ylabel("Iterations");
figure; plot(hset, recordingsc, "m");
xlabel("Resolution h"); ylabel("Calculations per iteration");
% % Plotting
% [X,Y] = meshgrid(x,y);
% figure; surf(X,Y,U); % 3D
% xlabel('x'); ylabel('y'); zlabel('U(x,y)');
% title("Resolution h = 0.01");
% %view(0,90);
% %saveas(gcf,'e1flat.png');
% view(-70,20);
% %saveas(gcf,'h001.png');
% 
% % figure; contour(X,Y,U); % from above
% % xlabel('x'); ylabel('y'); zlabel('U(x,y)');
% figure; meshc(X,Y,U) %mesh: contour + surf
% xlabel('x'); ylabel('y'); zlabel('U(x,y)');
% title("Resolution h = 0.01");
% view(-70,20);
% %saveas(gcf,'h001_mesh.png');