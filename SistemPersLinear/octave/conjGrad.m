function [x,numIter] = conjGrad(func,x,b,epsilon)
% Solves Ax = b by conjugate gradient method.
% USAGE: [x,numIter] = conjGrad(func,x,b,epsilon)
% INPUT:
% func = handle of function that returns the vector A*v
% x = starting solution vector
% b = constant vector in A*x = b
% epsilon = error tolerance (default = 1.0e-9)
% OUTPUT:
% x = solution vector
% numIter = number of iterations carried out
if nargin == 3
  epsilon = 1.0e-9;
end

n = length(b);

r = b - feval(func,x);
s = r;

for numIter = 1:n
  u = feval(func,s);
  alpha = dot(s,r)/dot(s,u);
  x = x + alpha*s;
  r = b - feval(func,x);
  if sqrt(dot(r,r)) < epsilon
    return
  else
    beta = -dot(r,u)/dot(s,u);
    s = r + beta*s;
  end
end
error('Too many iterations')
