function [xSol,ySol] = taylor(deriv,x,y,xStop,h)
% 4th-order Taylor series method of integration.
% USAGE: [xSol,ySol] = taylor(deriv,x,y,xStop,h)
% INPUT:
% deriv = handle of function that returns the matrix
% d = [dy/dx d^2y/dx^2 d^3y/dx^3 d^4y/dx^4].
% x,y = initial values; y must be a row vector.
% xStop = terminal value of x
% h = increment of x used in integration (h > 0).
% OUTPUT:
% xSol = x-values at which solution is computed.
% ySol = values of y corresponding to the x-values.

if size(y,1) > 1
  y = y'
end

% y must be row vector
xSol = zeros(2,1);
ySol = zeros(2,length(y));
xSol(1) = x;
ySol(1,:) = y;
k = 1;
while x < xStop
  h = min(h,xStop - x);
  d = feval(deriv,x,y);
  % Derivatives of [y]
  hh = 1;
  for j = 1:4
    hh = hh*h/j;
    % Build Taylor series
    % hh = h^j/j!
    y = y + d(j,:)*hh;
  end
  x = x + h; k = k + 1;
  xSol(k) = x; ySol(k,:) = y; % Store current soln.
end
