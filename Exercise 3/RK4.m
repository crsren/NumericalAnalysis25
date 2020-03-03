function [Vout, t, qc, qc_grad] = RK4(Vout, F2, F1, h, t, qc, qc_grad, N)
    for i=1:N 
        k1 = F1(t(i), qc(i), qc_grad(i));
        m1 = F2(t(i), qc(i), qc_grad(i));
        k2 = F1(t(i) + 0.5*h, qc(i) + (0.5*h*k1), qc_grad(i)+(0.5*h*m1));
        m2 = F2(t(i) + 0.5*h, qc(i) + (0.5*h*k1), qc_grad(i)+(0.5*h*m1));
        k3 = F1(t(i) + 0.5*h, qc(i) + (0.5*h*k2), qc_grad(i)+(0.5*h*m2));
        m3 = F2(t(i) + 0.5*h, qc(i) + (0.5*h*k2), qc_grad(i)+(0.5*h*m2));
        k4 = F1(t(i) + h, qc(i) + (k3*h), qc_grad(i) + (m3*h));
        m4 = F2(t(i) + h, qc(i) + (k3*h), qc_grad(i) + (m3*h));
        qc(i+1) = qc(i) + ((k1+k4)/6 + (k2+k3)/3)*h;
        qc_grad(i+1) = qc_grad(i) + ((m1+m4)/6 + (m2+m3)/3)*h;
        t(i+1) = t(i) + h;
        Vout(i+1) = qc_grad(i+1)*250;
    end
end 