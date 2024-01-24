% Using goldSearch to find x that minimizes:
% f(x) = 1.6*x^3 + 3*x^2 - 2*x
% s.t. constraint x >= 0

mu = 1.0;
func = @(x) 1.6*x^3 + 3.0*x^2 - 2.0*x + mu*max(0.0,-x)^2;

x = 1.0; % initial guess
[a,b] = goldBracket(func, x)
[xMin,fMin] = goldSearch(func, a, b)
