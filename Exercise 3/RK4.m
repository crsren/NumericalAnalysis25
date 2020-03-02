function [t, qc, qc_grad] = RK4(F, h, t, qc, qc_grad)
    for i=1: (length(t)-1)
        
        k1 = h*F(t(i), qc(i), qc_grad(i));
        k2 = h*F(t(i), qc(i) + h/2, qc_grad(i) + k1/2);
        k3 = h*F(t(i), qc(i) + h/2, qc_grad(i) + k2/2);
        k4 = h*F(t(i), qc(i) + h , qc_grad(i) + k3);
        
        
        qc_grad(i+1) = qc_grad(i) + (k1 + 2*k2 + 2*k3 + k4)/6;
        qc(i+1) = qc(i) + h*qc_grad(i);
        t(i+1) = t(i) + h;
    end
end