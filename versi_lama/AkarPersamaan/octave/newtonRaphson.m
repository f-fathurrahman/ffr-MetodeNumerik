function root = newtonRaphson( func, dfunc, a, b, tol )
% Newton-Raphson method combined with bisection

if nargin < 5
  tol = 1e6*eps
end

fa = feval(func, a);
if fa == 0
  root = a;
  return
end

fb = feval(func, b);
if fb == 0
  root = b;
  return
end

if fa*fb > 0.0
  error('Root is not bracketed in (a,b)')
end

x = (a + b)/2.0;

MaxIter = 30
for i = 1:MaxIter
  fx = feval(func, x);
  if abs(fx) < tol
    root = x;
    return
  end
  % Tighten brackets on the root
  if fa*fx < 0.0
    b = x;
  else
    a = x;
  end
  % Newton-Raphson step
  dfx = feval(dfunc, x);
  if abs(dfx) == 0
    dx = b - a;
  else
    dx = -fx/dfx;
  end
  %
  x = x + dx;
  %
  if (b-x)*(x-a) < 0.0
    dx = 0.5*(b-a);
    x = a + dx;
  end
  % convergence check
  if abs(dx) < tol*max(b,1.0)
    root = x;
    return
  end
end
root = NaN
