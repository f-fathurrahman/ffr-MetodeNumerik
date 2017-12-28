function [root,numIter] = newtonRaphson_simple(func, dfunc, x, tol)

if nargin < 5
  tol = 1.e6*eps;
end

MaxIter = 30;
for i = 1:MaxIter
  dx = -feval(func,x)/feval(dfunc,x);
  x = x + dx;
  if abs(dx) < tol
    root = x;
    numIter = i;
    return
  end
end
root = NaN