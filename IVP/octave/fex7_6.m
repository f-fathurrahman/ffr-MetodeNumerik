function F = fex7_6(x,y)
F = zeros(1,4);
F(1) = y(2);
F(2) = y(1)*y(4)^2 - 3.9860e14/y(1)^2;
F(3) = y(4);
F(4) = -2*y(2)*y(4)/y(1);