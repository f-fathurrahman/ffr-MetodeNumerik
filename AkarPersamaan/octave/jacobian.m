function [jac,f0] = jacobian(func,x)
% Returns the Jacobian matrix and f(x).
h = 1.0e-4;
n = length(x);
jac = zeros(n);
f0 = feval(func,x);
for i =1:n
  temp = x(i);
  x(i) = temp + h;
  f1 = feval(func,x);
  x(i) = temp;
  jac(:,i) = (f1 - f0)/h;
end
