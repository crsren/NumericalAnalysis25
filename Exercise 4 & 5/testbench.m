% SOR testing

% alpha < 1 : under relaxation
% alpha == 1: normal relaxation
% alpha > 1: over relaxation

close all ; clc

alpha = (0.8:0.1:2.1);
outputs = zeros(1, length(alpha));

for q = 1:length(alpha)
    disp("Testing alpha " + alpha(q));
    
    output = SOR(alpha(q));
    if(output < timeout)
        outputs(q) = output;
    else
        outputs(q) = max(outputs);
        break;
    end
end

figure; plot(alpha, outputs);
