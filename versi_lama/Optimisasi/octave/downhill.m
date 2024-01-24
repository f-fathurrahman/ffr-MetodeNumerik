function [xMin,fMin,k] = downhill(func,xStart,b,tol)
% Downhill simplex method for minimizing f(x1,x2,...,xn).
% USAGE: [xMin,nCycl] = downhill(func,xStart,b,tol)
% INPUT:
% func = handle of function to be minimized.
% xStart = coordinates of the starting point.
% b = initial edge length of the simplex (default = 0.1).
% tol = error tolerance (default = 1.0e-6).
% OUTPUT:
% xMin = coordinates of the minimum point.
% fMin = minimum value of f
% nCycl = number of cycles to convergence.

if nargin < 4
  tol = 1.0e-6
end

if nargin < 3
  b = 0.1
end

% Number of coords (design variables)
n = length(xStart);

% Coordinates of the n+1 vertices
x = zeros(n+1,n);

f = zeros(n+1,1);

% Values of 'func' at the vertices
% Generate starting simplex of edge length 'b'
x(1,:) = xStart;

for i = 2:n+1
  x(i,:) = xStart;
  x(i,i-1) = xStart(i-1) + b;
end

% Compute func at vertices of the simplex
for i = 1:n+1
  f(i) = func(x(i,:));
end

MaxIter = 500;
% Main loop
for k = 1:MaxIter
  % Find highest and lowest vertices
  [fLo,iLo]= min(f);
  [fHi,iHi] = max(f);
  % Compute the move vector d
  d = -(n+1)*x(iHi,:);
  for i = 1:n+1
    d = d + x(i,:);
  end
  d = d/n;
  % Check for convergence
  if sqrt(dot(d,d)/n) < tol
    xMin = x(iLo,:);
    fMin = fLo;
    return
  end
  % Try reflection
  xNew = x(iHi,:) + 2.0*d;
  fNew = func(xNew);
  if fNew <= f(iLo)
    % Accept reflection
    x(iHi,:) = xNew;
    f(iHi) = fNew;
    % Try expansion
    xNew = x(iHi,:) + d;
    fNew = func(xNew);
    if fNew <= f(iLo)
      % Accept expansion
      x(iHi,:) = xNew;
      f(iHi) = fNew;
      continue
    end
  else
    if fNew <= f(iHi)
      x(iHi,:) = xNew;
      f(iHi) = fNew;
      continue
    else
      % Try contraction
      xNew = x(iHi,:) + 0.5*d;
      fNew = func(xNew);
      if fNew <= f(iHi)
        % Accept contraction
        x(iHi,:) = xNew; f(iHi) = fNew;
        continue
      else
      % Shrinkage
        for i = 1:n+1
          if i ~= iLo
            x(i,:) = 0.5*(x(i,:) + x(iLo,:));
            f(i) = func(x(i,:));
          end
        end
      end
    end
  end
end

% End of main loop
error('Too many cycles in downhill')

