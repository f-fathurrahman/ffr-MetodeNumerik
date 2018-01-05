function [xSol,ySol] = bulStoer(dEqs,x,y,xStop,H,tol)
% Simplified Bulirsch-Stoer method for integration of y' = F(x,y).
% USAGE: [xSol,ySol] = bulStoer(dEqs,x,y,xStop,H,tol)
% INPUT:
% dEqs = handle of function that returns the first-order
%        differential equations F(x,y) = [dy1/dx,dy2/dx,...].
% x, y = initial values; y must be a row vector.
% xStop = terminal value of x.
% H = increment of x at which solution is stored.
% tol = per-step error tolerance (default = 1.0e-6).
% OUTPUT:
% xSol, ySol = solution at increments H.

if size(y,1) > 1
  y = y';
end

% y must be row vector
if nargin < 6
  tol = 1.0e-6;
end

n = length(y);

xSol = zeros(2,1);
ySol = zeros(2,n);
xSol(1) = x;
ySol(1,:) = y;

k = 1;
while x < xStop
  k = k + 1;
  H = min(H,xStop - x);
  y = midpoint(dEqs,x,y,x + H,tol);
  x = x + H;
  xSol(k) = x;
  ySol(k,:) = y;
end

