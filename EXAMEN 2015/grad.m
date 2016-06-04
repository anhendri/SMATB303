% =========================================================================
% Fonction pour evaluer le gradient de la fonction objectif en 1 point
% AS CRELOT, Juin 2015
%
% Input
%       pt : point de dimension n où evaluer la fonction
%            1 vecteur colonne  (nx1)
%       choix : indice indiquant la fonction a evaluer
%               1 pour la fonction quadratique
%               2 pour la fonction rosenbrock
%
% Output gsol: gradient evalue en pt
%              vecteur colonne (nx1)
% =========================================================================

function gsol = grad(pt,choix)
x = pt(1,1);
y = pt(2,1);
if (choix==1)
    gsol = [6*x+2*y-8; 2*x+6*y-8];
elseif (choix==2)
    gsol= [2 * (x-1) + 400* (x^2 -y) * (x);  -200 *(x^2 -y)  ];
end
end