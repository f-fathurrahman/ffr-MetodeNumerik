function [F] = F(W)

% ***** F REPRESENTS THE FUNCTION ON THE RIGHT-HAND SIDE OF EQUATION 
% (1.2.3). W STANDS FOR VELOCITY AND R FOR REYNOLDS NUMBER *****

global A B C D NU

R = abs(W)*D/NU;

% ***** STOKES FORMULA CANNOT BE USED WHEN R=0. IN THIS CASE THE VALUE
% ZERO IS ASSIGNED TO CD *****

if (R==0.)
    CD=0.;
end
if (R>0. && R<=1.)
    CD=24./R;
end
if (R>1. && R<=400.)
    CD=24./R^0.646;
end
if (R>400. && R<=3.0E+5)
    CD=0.5;
end
if (R>3.0E+5 && R<=2.0E+6)
    CD=3.66E-4*R^0.4275;
end
if (R>2.0E+6)
    CD=0.18;
end
F = (B - C*W*abs(W)*CD)/A;