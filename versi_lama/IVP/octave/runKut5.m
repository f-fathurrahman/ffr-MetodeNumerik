function [xSol,ySol] = runKut5(dEqs,x,y,xStop,h,eTol)
% 5th-order Runge--Kutta integration.
% USAGE: [xSol,ySol] = runKut5(dEqs,x,y,xStop,h,eTol)
% INPUT:
% dEqs = handle of function that specifies the
% 1st-order differential equations
% F(x,y) = [dy1/dx dy2/dx dy2/dx ...].
% x,y = initial values; y must be row vector.
% xStop = terminal value of x.
% h = trial value of increment of x.
% eTol = per-step error tolerance (default = 1.0e-6).
% OUTPUT:
% xSol = x-values at which solution is computed.
% ySol = values of y corresponding to the x-values.

if size(y,1) > 1
  y = y';
end

% y must be row vector
if nargin < 6
  eTol = 1.0e-6
end

n = length(y);

A = [0 1/5 3/10 3/5 1 7/8];
B = [0    0     0   0   0
     1/5  0     0   0   0
    3/40 9/40   0   0   0
    3/10 -9/10 6/5  0   0 
    -11/54  5/2 -70/27 35/27 0
    1631/55296 175/512 575/13824 44275/110592 253/4096];
C = [37/378 0 250/621 125/594 0 512/1771];
D = [2825/27648 0 18575/48384 13525/55296 277/14336 1/4];

% Initialize solution
xSol = zeros(2,1);
ySol = zeros(2,n);

xSol(1) = x;
ySol(1,:) = y;

stopper = 0;
k = 1;

for p = 2:5000
  % Compute K's from Eq. (7.18)
  K = zeros(6,n);
  K(1,:) = h*feval(dEqs,x,y);
  for i = 2:6
    BK = zeros(1,n);
    for j = 1:i-1
      BK = BK + B(i,j)*K(j,:);
    end
    K(i,:) = h*feval(dEqs, x + A(i)*h, y + BK);
  end
  % Compute change in y and per-step error from Eqs.(7.19) & (7.20)
  dy = zeros(1,n);
  E = zeros(1,n);
  for i = 1:6
    dy = dy + C(i)*K(i,:);
    E = E + (C(i) - D(i))*K(i,:);
  end
  e = sqrt(sum(E.*E)/n);
  % If error within tolerance, accept results and check for termination
  if e <= eTol
    y = y + dy;
    x = x + h;
    k = k + 1;
    xSol(k) = x;
    ySol(k,:) = y;
    if stopper == 1
      break
    end
  end
  % Size of next integration step from Eq. (7.24)
  hNext = 0.9*h*(eTol/e)^0.2;
  % Check if next step is the last one (works with positive and negative h)
  if (h > 0) == (x + hNext >= xStop )
    hNext = xStop - x; stopper = 1;
  end
  h = hNext;
end
