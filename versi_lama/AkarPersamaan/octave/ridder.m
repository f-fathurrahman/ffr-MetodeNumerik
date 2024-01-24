function root = ridder(func, x1, x2, tol)

if nargin < 4
  tol = 1.e6*eps;
end

f1 = func(x1);
if f1 == 0
  root = x1;
  return
end

f2 = func(x2);
if f2 == 0
  root = x2;
  return
end

if f1*f2 > 0
  error('Root is not bracketed in (a,b)')
end

for i = 0:30
  x3 = 0.5*(x1 + x2);
  f3 = func(x3);
  if f3 == 0.0
    root = x3;
    return
  end
  % Ridder formula
  s = sqrt(f3^2 - f1*f2);
  if s == 0
    root = NaN;
    return
  end
  %
  dx = (x3 - x1)*f3/s;
  if (f1 - f2) < 0
    dx = -dx;
  end
  x4 = x3 + dx;
  f4 = func(x4);
  %
  if i > 0
    if abs(x4 - xOld) < tol*max(abs(x4),1.0)
      root = x4;
      return
    end
  end
  xOld = x4;
  %
  if f3*f4 > 0
    if f1*f4 < 0
      x2 = x4;
      f2 = f4;
    else
      x1 = x4;
      f1 = f4;
    end
  else
    x1 = x3;
    x2 = x4;
    f1 = f3;
    f2 = f4;
  end
end
root = NaN;