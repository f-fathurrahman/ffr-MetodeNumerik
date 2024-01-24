func = @(t) (2* t^2 * cos(t^2));

N = 10;
Integral = gaussQuad( func, 0, sqrt(pi), N );

fprintf('Integral = %18.10f\n', Integral)
fprintf('N        = %d\n', N)
