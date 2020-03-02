function [t, qc, qc_grad] = RK4(F, h, t, qc, qc_grad)
    for i=1: (length(t)-1)
        
        k1 = qc_grad(i);
        m1 = F(t(i), qc(i), qc_grad(i));
        
        k2 = qc_grad(i) + (m1*h/2);
        m2 = F(t(i) + 0.5*h, qc(i) + (0.5*h*k1), qc_grad(i)+(0.5*h*m1));
        
        k3 = qc_grad(i) + (m2*h);
        m3 = F(t(i) + 0.5*h, qc(i) + (0.5*h*k2), qc_grad(i)+(0.5*h*m2));
        
        k4 = qc_grad(i) + (m3);
        m4 = F(t(i) + h, qc(i) + (k3*h), qc_grad(i) + (m3*h));
        
        qc_grad(i+1) = qc_grad(i) + (k1+2*k2+2*k3+k4)*h/6;
        qc(i+1) = qc(i) + (m1+2*m2+2*m3+m4)*h/6;
        t(i+1) = t(i) + h;
    end
end