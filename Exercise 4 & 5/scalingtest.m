% Testing how iterations change given different constant boundary
% conditions

% Set boundary conditions
Ui_min = @(x) x;%(x-100)^2;
Ui_max = @(x) x;%(x-100)^2;
Umin_j = @(y) y;%(y-100)^2;
Umax_j = @(y) y;%(y-100)^2;

height = 1:10:200;
itrec = zeros(1, length(height));

for index=1:length(height)
    
    h = 0.01; n = 1; m = 1;
    x = 0:h:n;
    y = 0:h:m;
    U = zeros(length(x),length(y));

    for i = 1 : length(x)
        U(i ,1) = height(index)*Ui_min(x(i));
        U(i,length(y)) = height(index)*Ui_max(x(i));
    end
    for j = 1 : length(y)
        U(1, j) = height(index)*Umin_j(y(j));
        U(length(x), j) = height(index)*Umax_j(y(j));
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
    
    disp(height(index) + ": Done after " + iterations + " iterations");

    itrec(index) = iterations;
end

figure; plot(height, itrec, "r*-");
xlabel("BC height"); ylabel("Iterations");
title("Effect of different constant BC");
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