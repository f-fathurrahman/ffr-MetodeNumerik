function [xMin, fMin] = goldSearch(func, a, b, tol)

if nargin < 4
  tol = 1.0e-6;
end

nIterMax = ceil(-2.078087*log(tol/abs(b-a)));
R = 0.618033989;
C = 1.0 - R;

% First telescoping
x1 = R*a + C*b;
x2 = C*a + R*b;

f1 = feval(func, x1);
f2 = feval(func, x2);

% Main loop
for i = 1:nIterMax
  if f1 > f2
    a = x1;
    x1 = x2;
    f1 = f2;
    x2 = C*a + R*b;
    f2 = feval(func, x2);
  else
    b = x2;
    x2 = x1;
    f2 = f1;
    x1 = R*a + C*b;
    f1 = feval(func, x1);
  end
end

% Determine the output
if f1 < f2
  fMin = f1;
  xMin = x1;
else
  fMin = f2;
  xMin = x2;
end
