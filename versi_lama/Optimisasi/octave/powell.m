function [xMin,fMin,nCyc] = powell(func,x,h,tol)
% Powell’s method for minimizing f(x1,x2,...,xn).
% USAGE: [xMin,fMin,nCyc] = powell(h,tol)
% INPUT:
% func = handle of function that returns f.
% x = starting point
% h = initial search increment (default = 0.1).
% tol = error tolerance (default = 1.0e-6).
% OUTPUT:
% xMin = minimum point
% fMin = minimum value of f
% nCyc = number of cycles to convergence

if nargin < 4
  tol = 1.0e-6;
end

if nargin < 3
  h = 0.1;
end

% x must be column vector
if size(x,2) > 1
  x = x';
end

% Number of design variables
n = length(x);

df = zeros(n,1); % Decreases of f stored here

u = eye(n); % Columns of u store search directions v

MaxIter = 30;
for j = 1:MaxIter
  xOld = x;
  fOld = feval(func,xOld);
  % First n line searches record the decrease of f
  for i = 1:n
    v = u(1:n,i);
    fLine = @(s) fLine_full(func, x, v, s);
    [a,b] = goldBracket(fLine,0.0,h);
    [s,fMin] = goldSearch(fLine,a,b);
    df(i) = fOld - fMin;
    fOld = fMin;
    x = x + s*v;
  end
  % Last line search in the cycle
  v = x - xOld;
  fLine = @(s) fLine_full(func, x, v, s);
  [a,b] = goldBracket(fLine,0.0,h);
  [s,fMin] = goldSearch(fLine,a,b);
  x = x + s*v;
  % Check for convergence
  if sqrt(dot(x-xOld,x-xOld)/n) < tol
    xMin = x;
    nCyc = j;
    return
  end
  % Identify biggest decrease of f & update search directions
  iMax = 1; dfMax = df(1);
  for i = 2:n
    if df(i) > dfMax
      iMax = i; dfMax = df(i);
    end
  end
  
  for i = iMax:n-1
    u(1:n,i) = u(1:n,i+1);
  end
  u(1:n,n) = v;
end

error('Powell method did not converge')

endfunction


function z = fLine_full(func,x,v,s) % F in the search direction v
z = feval(func, x + s*v);
endfunction

