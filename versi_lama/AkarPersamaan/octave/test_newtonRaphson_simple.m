func = @(x) (x^3 - 10*x^2 + 5);
dfunc = @(x) (3*x^2 - 20*x);

[x0, numIter] = newtonRaphson_simple( func, dfunc, 0.7 )


% Example 4.7 (Newton--Raphson method)
func = @(x) (x^4 - 6.4*x^3 + 6.45*x^2 + 20.538*x - 31.752);
dfunc = @(x) (4*x^3 - 19.2*x^2 + 12.9*x + 20.538);

xStart = 2.0;
[x0,numIter] = newtonRaphson_simple(func, dfunc, xStart)