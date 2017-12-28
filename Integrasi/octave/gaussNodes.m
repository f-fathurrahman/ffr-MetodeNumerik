function [x,A] = gaussNodes(N,tol,maxIter)

if nargin < 3
  maxIter = 50;
end

if nargin < 2
  tol = 1.e4*eps;
end

A = zeros(N,1);
x = zeros(N,1);


% Truncate fractional portion of X and return the integer portion.
% This is equivalent to rounding towards zero.
nRoots = fix(N+1)/2;

%fprintf('nRoots = %d\n', nRoots)

for i = 1:nRoots
  t = cos(pi*(i-0.25)/(N+0.5)); % approximate roots
  % Use Newton's root-finding method
  for j = i:maxIter
    [p, dp] = my_legendre(t,N);
    dt = -p/dp;
    t = t + dt;
    if abs(dt) < tol
      x(i) = t;
      x(N-i+1) = -t;
      A(i) = 2/(1-t^2)/dp^2;
      A(N-i+1) = A(i);
      break
    end
  end
end

