function [x,numIter,omega] = gaussSeidel(func,x,maxIter,epsilon)
% Solves Ax = b by Gauss-Seidel method with relaxation.
% USAGE: [x,numIter,omega] = gaussSeidel(func,x,maxIter,epsilon)
% INPUT:
% func = handle of function that returns improved x using
%        the iterative formula
% x = starting solution vector
% maxIter = allowable number of iterations (default is 500)
% epsilon = error tolerance (default is 1.0e-9)
% OUTPUT:
% x = solution vector
% numIter = number of iterations carried out
% omega = computed relaxation factor
if nargin < 4
  epsilon = 1.0e-9
end

if nargin < 3
  maxIter = 500
end

k = 10;
p = 1;
omega = 1;

for numIter = 1:maxIter
  xOld = x;
  x = feval(func,x,omega);
  dx = sqrt(dot(x - xOld,x - xOld));
  if dx < epsilon
    return
  end
  if numIter == k
    dx1 = dx;
  end
  if numIter == k + p
    omega = 2/(1 + sqrt(1 - (dx/dx1)^(1/p)));
  end
end
error('Too many iterations')
