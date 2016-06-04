% =========================================================================
% Fonction pour evaluer la fonction objectif en 1 ou plusieurs points
% AS CRELOT, Juin 2015
%
% Input
%       pt : point(s) de dimension n où evaluer la fonction
%            1 vecteur colonne si un seul point (nx1)
%            1 matrice (nxm) si m points (chaque point est stocke sur une
%            colonne)
%       choix : indice indiquant la fonction a evaluer
%               1 pour la fonction quadratique
%               2 pour la fonction rosenbrock
%
% Output fsol: fonction evaluee en pt
%              scalaire si un seul point
%              vecteur ligne (dimension 1xm) si m points
% =========================================================================

function fsol = fct(pt,choix)
x=pt(1,:);
y=pt(2,:);
if (choix==1)
    fsol = 2.*(x+y-2).^2 + (x-y).^2;
elseif (choix==2)
    fsol= 100*(y-x.*x).^2+(1-x).^2;
end
end