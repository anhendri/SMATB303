% =========================================================================
% Algorithme Recherche linéaire pour la longueur de pas
% Conditions de Wolfe
% Algo 3.5 p.60 Nocedal and Wright
% AS CRELOT, Juin 2015
%
% Input
%       d : direction de descente
%       x0 : itere courrant (a partir duquel on realise la recherche)
%       alphamax : pas maximum
%       choix : indice indiquant la fonction objectif utilisee
%               1 pour la fonction quadratique
%               2 pour la fonction rosenbrock
%
% Output alphas : pas verifiant les conditions de wolfe
% =========================================================================

function alphas = linesearchWolfe(d,x0,alphamax,choix)

c1 = 1e-4;
c2 = 0.5;
alpha0 = 0;
alphap = alpha0;         % alphap is alpha_{i-1}
alphax = alphamax/2;     % alphax is alpha_i
i=1;

fx0 = fct(x0,choix);     %phi(0)
gx0 = grad(x0,choix);
 
fxp = fx0;               %phi(alpha_{i-1})

while (true)
  xx = x0 + alphax*d;
  fxx = fct(xx,choix);   %phi(alpha_i)
  gxx = grad(xx,choix);
  
  if (fxx > fx0 + c1*alphax*gx0'*d) || ((i > 1) && (fxx >= fxp)),
    alphas = zoom(x0,d,alphap,alphax,choix);
    return;
  end
  if abs(gxx'*d) <= -c2*gx0'*d
    alphas = alphax;
    return;
  end
  if gxx'*d >= 0,
    alphas = zoom(x0,d,alphax,alphap,choix);
    return;
  end
  
  alphap = alphax;
  fxp = fxx;
  alphax = min(alphamax, alphax*3);
  i = i+1;
end

