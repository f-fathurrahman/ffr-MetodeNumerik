function F = fex7_11(x,y)
F = zeros(1,2);
R = 1.0;
L = 2.0;
C = 0.45;
% y(1) = q
% y(2) = i
E = 9.0;
F(1) = y(2);
F(2) = ( -R*y(2) - y(1)/C + E )/L;