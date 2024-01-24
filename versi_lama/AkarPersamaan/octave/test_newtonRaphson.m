func = @(x) (x^3 - 10*x^2 + 5);
dfunc = @(x) (3*x^2 - 20*x);

x0 = newtonRaphson( func, dfunc, 0.6, 0.8 )