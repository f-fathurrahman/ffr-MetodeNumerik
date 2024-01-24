function F = fex7_8(x,y)
F = zeros(1,2);
F(1) = y(2);
F(2) = -9.80665 + 65.351e-3 * y(2)^2 * exp(-10.53e-5 * y(1));