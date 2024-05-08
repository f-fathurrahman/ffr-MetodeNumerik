function [RES] = F(T,VEC);

% ***** IT REPRESENTS THE RIGHT-HAND SIDE OF EQUATION (1.3.12) *****

global ALPHA0 BETA PI U

Z = VEC(1);
V = VEC(2);

F = -(2.*PI)^2*Z + BETA*U*sqrt(U*U+V*V)*CL(V);

RES = [V;F];
