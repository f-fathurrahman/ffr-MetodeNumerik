function y = midpoint(dEqs,x,y,xStop,tol)
% Modified midpoint method for integration of y' = F(x,y).
% USAGE: y = midpoint(dEqs,xStart,yStart,xStop,tol)
% INPUT:
% dEqs = handle of function that returns the first-order
%        differential equations F(x,y) = [dy1/dx,dy2/dx,...].
% x, y = initial values; y must be a row vector.
% xStop = terminal value of x.
% tol = per-step error tolerance (default = 1.0e-6).
% OUTPUT:
% y = y(xStop).

  if size(y,1) > 1
    y = y';
  end

  % y must be row vector
  if nargin <5
    tol = 1.0e-6;
  end

  kMax = 21;

  n = length(y);
  r = zeros(kMax,n);

  % Storage for Richardson extrapolation.
  % Start with two integration steps.
  nSteps = 2;
  r(1,1:n) = mid(dEqs,x,y,xStop,nSteps);
  rOld = r(1,1:n);
  for k = 2:kMax
    % Increase the number of steps by 2 & refine results by Richardson extrapolation.
    nSteps = 2*k;
    r(k,1:n) = mid(dEqs,x,y,xStop,nSteps);
    r = richardson(r,k,n);
    % Check for convergence.
    dr = r(1,1:n) - rOld;
    e = sqrt(dot(dr,dr)/n);
    if e < tol
      y = r(1,1:n);
      return;
    end
    rOld = r(1,1:n);
  end
  error('Midpoint method did not converge')
endfunction

function r = richardson(r,k,n)
  % Richardson extrapolation.
  for j = k-1:-1:1
    c =(k/(k-1))^(2*(k-j));
    r(j,1:n) =(c*r(j+1,1:n) - r(j,1:n))/(c - 1.0);
  end
  return
endfunction


function y = mid(dEqs,x,y,xStop,nSteps)
  % Midpoint formulas.
  h = (xStop - x)/nSteps;
  y0 = y;
  y1 = y0 + h*feval(dEqs,x,y0);
  for i = 1:nSteps-1
    x = x + h;
    y2 = y0 + 2.0*h*feval(dEqs,x,y1);
    y0 = y1;
    y1 = y2;
  end
  y = 0.5*(y1 + y0 + h*feval(dEqs,x,y2));
endfunction