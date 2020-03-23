% SOR testing

% alpha < 1 : under relaxation
% alpha == 1: normal relaxation
% alpha > 1: over relaxation

%close all ; clc

% Initializing vectors to test different α
alpha = horzcat((0.8:0.1:1.4),(1.42:0.02:2));
outputs = zeros(1, length(alpha));
timeout = 10000;

% d2u/dx2 + d2u/dy2 = g(x,y)
g = @(x,y) 0;

% Define grid
h = 0.01; n = 1; m = 1;

% Set boundary conditions: 1,0; 10,0; x; 0.5*(2*x-1)^2; 0.5*(2*x-1)^4; sin(2*pi*x);
Ui_min = @(x) x;
Ui_max = @(x) x;
Umin_j = @(y) y;
Umax_j = @(y) y;

% Run
run = true;
for q = 1:(length(alpha)-1)
    if run

        disp("Testing alpha " + alpha(q));
        output = SOR(alpha(q),h,n,m,Ui_min,Ui_max,Umin_j,Umax_j,g);
        
        if(output > timeout)
            plotlim = round(max(outputs/100))*100;
            disp("gude morge");
            run = false;
        end
        outputs(q) = output;
    end
end

% Adding alpha = 2 (does not terminate)
outputs(length(alpha)) = 10000;

%% Plotting iterations vs relaxation factor
hold on;
plot(alpha, outputs); ylim([0,3000]); xlim([0.8,2.01]);
xlabel("Relaxation factor α");
ylabel("Iterations");
title("SOR with varying α for different BC");


function [yy] = discon(xx)
    if(xx < 0.25)
        yy = 0;
    elseif(xx < 0.5) 
        yy = 1;
    elseif(xx < 0.75)
        yy = 0;
    else
        yy = 1;
    end
end

%% Successive over-relaxation (SOR.m)

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
    
    % Apply SOR method on all interior points until the residue r < Ɛ
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
    
    %% Plotting solutions
    %[X,Y] = meshgrid(x,y);
    %figure; surf(X,Y,U); % 3D
    %figure; contour(X,Y,U); % from above
    %figure; meshc(X,Y,U) %mesh: contour + surf
    %xlabel('x'); ylabel('y'); zlabel('U(x,y)'); title("α = " + a);
end
