function root = bisection(func, x1, x2, do_filter, tol)
% do_filter: use singularity filter or not

if nargin < 5
  tol = 1.0*eps;
end

if nargin < 4
  do_filter = false;
end

f1 = feval(func, x1);
if f1 == 0.0
  root = x1;
  return
end

f2 = feval(func, x2);
if f2 == 0.0
  root = x2;
  return
end

if f1*f2 > 0
  error('Root is not bracketed in (x1,x2)')
end

MaxIter = ceil( log(abs(x2-x1)/tol)/log(2.0) );
fprintf('MaxIter = %d\n', MaxIter)

for i = 1:MaxIter
  x3 = 0.5*(x1 + x2);
  f3 = feval(func, x3);
  %
  if( (do_filter == true) && (abs(f3) > abs(f1)) && (abs(f3) > abs(f2) ) )
    root = NaN;
    return
  end
  %
  if f3 == 0.0
    root = x3;
    return
  end
  %
  if f2*f3 < 0.0
    x1 = x3;
    f1 = f3;
  else
    x2 = x3;
    f2 = f3;
  end
end

root = 0.5*(x1 + x2);
