function [a,b] = goldBracket(func, x1, h)
% Brackets the minimum point of f(x)

if nargin < 3
  h = 0.1;
end

c = 1.618033989;

f1 = feval(func, x1);

x2 = x1 + h;
f2 = feval(func, x2);

% Determine downhill direction and change sign of h if needed
if f2 > f1
  h = -h;
  x2 = x1 + h;
  f2 = feval(func, x2);
  % check if minimum is between x1 - h and x1 + h
  if f2 > f1
    a = x2;
    b = x1 - h;
    return
  end
end

% search loop
MaxIter = 100;
for i = 1:100
  h = c*h;
  x3 = x2 + h;
  f3 = feval(func, x3);
  if f3 > f2
    a = x1;
    b = x3;
    return
  end
  x1 = x2;
  f1 = f2;
  x2 = x3;
  f2 = f3;
end
error('goldBracket did not find a minimum')

