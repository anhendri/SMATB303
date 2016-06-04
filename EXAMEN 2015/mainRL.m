function [] = mainRL()
    format long e
    dd = menu('Direction de descente','Plus forte pente','Newton',      ...
              'Quasi-Newton BFSG','Quasi-Newton SR1');

    lp = menu('Longueur de pas','Exacte Quadratique',                   ...
              'Recherche lineaire + Wolfe');

    cf = menu('Choix de la fonction a minimiser','Fonction quadratique',...
              'Rosenbrock');

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

    [niter,nitermax,tol, dk] = deal(1,10000,1E-10,[1,1]);                  % INITITATION DES VARIABLES
    B = eye(2);
    while(norm(dk)> tol && niter < nitermax)
        dk = -grad(x0,cf);                                                 % CALCUL DU GRADIENT EN Xk
        ddk = hessien(x0,cf);                                              % CALCUL DE LA MATRICE HESSIENNE EN Xk
        ak = (-dk)'*(-dk)/((-dk)'*(ddk*(-dk)));                            % CALCUL DE ALPHAk
        if(dd==2)
            dk = ddk\dk;
        elseif(dd==3)
            dk = B\dk;
            s = ak*dk;
            y = grad(x0+s,cf)-grad(x0,cf);
            if(s'*y>=0.2*s'*B*s)
                theta = 1;
            else
                theta = (0.8*s'*B*s)/(s'*B*s-s'*y);
            end
            t = theta*y+(1-theta)*B*s;
            B = B-(B*(s*s')*B)/(s'*B*s)+(t*t')/(s'*t);
        elseif(dd==4)
            dk = B\dk;
            s = ak*dk;
            y = grad(x0+s,cf)-grad(x0,cf);
            if(abs(s'*(y-B*s))<1E-8*norm(s)*norm(y-B*s))
                B = B + (y-B*s)*(y-B*s)'/((y-B*s)'*s);
            end
        end
     %   if(lp==2)
     %       ak = linesearchWolfe(dk,x0,100,cf); 
     %   end
        x0 = x0+ak*dk;                                                     % CALCUL DE Xk+1
        niter = niter + 1;                                                % INCREMENTATION DU NOMBRE D'ITERATIONS
    end
    disp(x0)
    niter
end