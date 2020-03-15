% SOR testing

% alpha < 1 : under relaxation
% alpha == 1: normal relaxation
% alpha > 1: over relaxation

%close all ; clc
figure;

alpha = horzcat((0.8:0.1:1.4),(1.42:0.02:2));
outputs = zeros(1, length(alpha));
timeout = 10000;

% d2u/dx2 + d2u/dy2 = g(x,y)
g = @(x,y) 0; % !! Poisson eqn other than Laplace not working yet !!

% Define grid
h = 0.01; n = 1; m = 1;

% Set boundary conditions: 1,0; 10,0; x; 0.5*(2*x-1)^2; 0.5*(2*x-1)^4; sin(2*pi*x);
Ui_min = @(x) x;
Ui_max = @(x) e^x;
Umin_j = @(y) e^y;
Umax_j = @(y) e^y;

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

%simulating alpha = 2
outputs(length(alpha)) = 10000;

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