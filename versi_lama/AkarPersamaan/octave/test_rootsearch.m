func = @(x) (x^3 - 10*x^2 + 5);

dx = 0.1;
a = 0;
b = 1;

[x1,x2] = rootsearch( func, a, b, dx )
