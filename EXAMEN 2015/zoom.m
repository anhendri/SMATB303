% =========================================================================
% Algorithme zoom utilise par la recherche lineaire avec conditions de
% Wolfe 
% Algo 3.6 p.61 Nocedal and Wright
% AS CRELOT, Juin 2015
%
% Fonction appelee par linsesearchWolfe
%
% Input
%       x0 : itere courrant (a partir duquel on realise la recherche)
%       d : direction de descente
%       alphal : borne inf pour le pas
%       alphah : borne sup pour le pas
%       choix : indice indiquant la fonction objectif utilisee
%               1 pour la fonction quadratique
%               2 pour la fonction rosenbrock
%
% Output alphas : pas verifiant les conditions de wolfe
% =========================================================================

function alphas = zoom(x0,d,alphal,alphah,choix)

c1 = 1e-4;
c2 = 0.5;

fx0 = fct(x0,choix);     %phi(0)
gx0 = grad(x0,choix);

k=0;

while (k<10000)
    k=k+1;
    
    alphax = (alphal+alphah)/2; %alpha_j
    
    xx = x0 + alphax*d;
    fxx = fct(xx,choix);   %phi(alpha_j)
    gxx = grad(xx,choix);
    
    xl = x0 + alphal*d;
    fxl = fct(xl,choix);   %phi(alphal)
    gxl = grad(xl,choix);
    
    if ((fxx > fx0 + c1*alphax*gx0'*d) || (fxx >= fxl))
        alphah = alphax;
    else
        if (abs(gxx'*d) <= -c2*gx0'*d)
            alphas = alphax;
            return;
        end
        if (gxx'*d*(alphah-alphal) >= 0)
            alphah = alphal;
        end
        alphal = alphax;
    end
end
alphas = alphax;



