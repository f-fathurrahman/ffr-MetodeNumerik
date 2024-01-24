function [xSol,ySol] = runKut4(dEqs,x,y,xStop,h)
% 4th-order Runge--Kutta integration.
% USAGE: [xSol,ySol] = runKut4(dEqs,x,y,xStop,h)
% INPUT:
% dEqs = handle of function that specifies the
% 1st-order differential equations
% F(x,y) = [dy1/dx dy2/dx dy2/dx ...].
% x,y = initial values; y must be row vector.
% xStop = terminal value of x.
% h = increment of x used in integration.
% OUTPUT:
% xSol = x-values at which solution is computed.
% ySol = values of y corresponding to the x-values.
if size(y,1) > 1
  y = y';
end

% y must be row vector
xSol = zeros(2,1);
ySol = zeros(2,length(y));
xSol(1) = x;
ySol(1,:) = y;

i = 1;
while x < xStop
  i = i + 1;
  h = min(h,xStop - x);
  K1 = h*feval(dEqs,x,y);
  K2 = h*feval(dEqs,x + h/2,y + K1/2);
  K3 = h*feval(dEqs,x + h/2,y + K2/2);
  K4 = h*feval(dEqs,x+h,y + K3);
  y = y + (K1 + 2*K2 + 2*K3 + K4)/6;
  x = x + h;
  % Store current solution
  xSol(i) = x;
  ySol(i,:) = y; 
end
