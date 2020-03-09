% SOR testing

% alpha < 1 : under relaxation
% alpha == 1: normal relaxation
% alpha > 1: over relaxation

close all ; clc

alpha = 1;%(0.5:0.1:2.1);
outputs = zeros(1, length(alpha));

% d2u/dx2 + d2u/dy2 = g(x,y)
g = @(x,y) 0; % !! Poisson eqn other than Laplace not working yet !!

% Define grid
h = 0.1; n = 100; m = 100;

% Set boundary conditions
Ui_min = @(x) x;
Ui_max = @(x) x;
Umin_j = @(y) y;
Umax_j = @(y) y;

for q = 1:length(alpha)
    disp("Testing alpha " + alpha(q));
    
    output = SOR(alpha(q),h,n,m,Ui_min,Ui_max,Umin_j,Umax_j,g);
    if(output < timeout)
        outputs(q) = output;
    else
        outputs(q) = max(outputs);
        break;
    end
end

figure; plot(alpha, outputs);