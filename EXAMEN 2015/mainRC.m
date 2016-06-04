function [] = mainRC()
    format long e
    clc

    disp('Direction de descente :');
    disp('-----------------------');
    disp('1 : Plus forte pente');
    disp('2 : Newton');
    disp('3 : Quasi-Newton SR1');
    disp('4 : Quasi-Newton BFGS');
    dd = input('Choix = ');

    disp('Choix de la fonction a minimiser :');
    disp('----------------------------------');
    disp('1 : Fonction quadratique');
    disp('2 : Fonction de Rosenbrock');
    cf = input('Choix = ');          

    if (cf==1)
        [x,y]=meshgrid(0:0.1:4,-2:0.1:2);
        Z=2.*(x+y-2).^2 + (x-y).^2;
        x0=[3;1.5];
    elseif (cf==2)
        [x,y]=meshgrid(-2:0.1:2,-2:0.1:4);
        Z=100.*(y-x.^2).^2+(1-x).^2;
        x0=[-1;1];
    else
        disp('erreur de choix');
        return;
    end
    figure;
    contour(x,y,Z,20);
    hold on;

    [niter,nitermax,tol, dk] = deal(1,10000,1E-10,[1,1]);                  % INITIATION DES VARIABLES
    if(dd==1)
        B = zeros(2);
    elseif(dd==2)
        B = hessien(x0,cf);
    else
        B = eye(2);
    end
    [eta1,eta2,alpha1,alpha2,delta] = deal(0.01,0.9,2.5,0.25,1);
    
    while(norm(dk)> tol && niter < nitermax)
        dk = grad(x0,cf);
        pk = SteighaugToint(dk,B,delta);
        rhok = -(fct(x0,cf)-fct(x0+pk,cf))/(dk'*pk + 1/2*pk'*B*pk);
        
        if(rhok>=eta1)
            x0 = x0 + pk;
            if(dd==1)
                B = zeros(2);
            elseif(dd==2)
                B = hessien(x0,cf);
            elseif(dd>=2)
                s = pk;
                y = grad(x0,cf)-grad(x0-pk,cf);
                if(abs(s'*(y-B*s))>=1E-8*norm(s)*norm(y-B*s))
                    B = B + (y-B*s)*(y-B*s)'/((y-B*s)'*s);
                end
            end
        end
        if(rhok >= eta2)
            delta = max([alpha1*norm(pk),delta]);
        elseif(rhok<eta1)
            delta = alpha2 * norm(pk);
        end
        niter = niter + 1;
    end
    disp(x0)
end