function [qc, qc_grad] = RK4(F1, F2, h, t, qc, qc_grad)
    for i=1: (length(t)-1) 
        k1 = F1(t(i), qc(i), qc_grad(i));
        m1 = F2(t(i), qc(i), qc_grad(i));
        k2 = F1(t(i) + 0.5*h, qc(i) + (0.5*h*k1), qc_grad(i)+(0.5*h*m1));
        m2 = F2(t(i) + 0.5*h, qc(i) + (0.5*h*k1), qc_grad(i)+(0.5*h*m1));
        k3 = F1(t(i) + 0.5*h, qc(i) + (0.5*h*k2), qc_grad(i)+(0.5*h*m2));
        m3 = F2(t(i) + 0.5*h, qc(i) + (0.5*h*k2), qc_grad(i)+(0.5*h*m2));
        k4 = F1(t(i) + h, qc(i) + (k3*h), qc_grad(i) + (m3*h));
        m4 = F2(t(i) + h, qc(i) + (k3*h), qc_grad(i) + (m3*h));
        qc(i+1) = qc(i) + (k1+2*k2+2*k3+k4)*h/6;
        qc_grad(i+1) = qc_grad(i) + (m1+2*m2+2*m3+m4)*h/6;
        t(i+1) = t(i) + h;
    end
end