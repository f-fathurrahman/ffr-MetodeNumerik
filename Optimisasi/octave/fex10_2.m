function S = fex10_2(y)
% Function used in Example 10.2
B = 48.0;
H = 60.0;

a = B*(H - y)/H;
b = (B - a)/2.0;

A = (B + a)*y/2.0;

Q = (a*y^2)/2.0 + (b*y^2)/3.0;
d = Q/A;
c = y - d;

I = (a*y^3)/3.0 + (b*y^3)/6.0;

Ibar = I - A*d^2; S = -Ibar/c;