function [CD] = CD(V)

% ***** CD IS THE APPROXIMATE DRAG COEFFICIENT OF A SPHERE *****

global A B C D G L NU

RE = abs(V)*D/NU;

if (RE==0.)
    CD=0.;
end
if (RE>0. && RE<=1.)
    CD=24./RE;
end
if (RE>1. && RE<=400.)
    CD=24./RE^0.646;
end
if (RE>400. && RE<=3.0E+5)
    CD=0.5;
end
if (RE>3.0E+5 && RE<=2.0E+6)
    CD=3.66E-4*RE^0.4275;
end
if (RE>2.0E+6)
    CD=0.18;
end
