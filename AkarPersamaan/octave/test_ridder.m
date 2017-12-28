func = @(x) (x^3 - 10*x^2 + 5);

x0 = ridder( func, 0.6, 0.8 )