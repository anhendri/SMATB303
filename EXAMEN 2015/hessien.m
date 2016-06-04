% =========================================================================
% Fonction pour evaluer le hessien de la fonction objectif en 1 point
% AS CRELOT, Juin 2015
%
% Input
%       pt : point de dimension n où evaluer la fonction
%            1 vecteur colonne  (nx1)
%       choix : indice indiquant la fonction a evaluer
%               1 pour la fonction quadratique
%               2 pour la fonction rosenbrock
%
% Output hsol: hessien evalue en pt
%              matrice (nxn)
% =========================================================================

function hsol = hessien(pt,choix)
x = pt(1,1);
y = pt(2,1);
if (choix==1)
    hsol = [6  2 ;2  6 ];
elseif (choix==2)
    hsol= [-400*y+1200*x.^2+2  -400*x ; -400*x 200  ];
end
end