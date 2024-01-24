function root = linInterp(func,x1,x2)
% Finds the zero of the linear function f(x) by straight
% line interpolation between x1 and x2.
% func = handle of function that returns f(x).
f1 = feval(func,x1);
f2 = feval(func,x2);
root = x2 - f2*(x2 - x1)/(f2 - f1);