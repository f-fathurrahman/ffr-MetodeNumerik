function [RES] = FN(T,VEC);

global A B C D G L NU

Y = VEC(1);
V = VEC(2);

FN = (-B*sin(Y/L) - C*V*abs(V)*CD(V)) / A;

RES = [V;FN];