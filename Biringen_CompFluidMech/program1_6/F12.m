function [RES] = F12(T,VEC);

global A B

X = VEC(1);
U = VEC(2);
Y = VEC(3);
V = VEC(4);

% ***** IT REPRESENTS THE RIGHT-HAND SIDE OF EQUATION (1.7.4) *****

F1 = -A * sqrt(U*U+V*V) * (B*U+V);

% ***** IT REPRESENTS THE RIGHT-HAND SIDE OF EQUATION (1.7.5) *****

F2 = -1. + A * sqrt(U*U+V*V) * (U-B*V);

RES = [U;F1;V;F2];
