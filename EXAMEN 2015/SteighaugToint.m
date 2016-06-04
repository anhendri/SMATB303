function pk = SteighaugToint(g,B,deltak)
    pk  = zeros(2,1);
    z   = pk;
    rk  = g;
    dk  = -rk;
    eta = 1e-6;

    while(norm(rk) >= eta)
        if((dk'*B*dk) <= 0)
            tau = (sqrt((z'*dk)^2+dk'*dk*(deltak^2-z'*z))-z'*dk)/(dk'*dk);
            pk  = z+tau*dk;
            return 
        end
        alpha = (rk'*rk)/(dk'*B*dk);
        z = z+alpha*dk;
        if(norm(z) >= deltak)
            tau = (sqrt((z'*dk)^2+dk'*dk*(deltak^2-z'*z))-z'*dk)/(dk'*dk);
            pk  = z+tau*dk;
            return 
        end
        rj = rk;
        rk = rk+alpha*B*dk;
        if(norm(rk)<eta)
            pk = z;
            return 
        end
        beta = (rk'*rk)/(rj'*rj);
        dk = beta*dk - rk;
    end
end