function [c,d,e] = LUdec3(c,d,e)
% LU decomposition of tridiagonal matrix A = [c\d\e].
% USAGE: [c,d,e] = LUdec3(c,d,e)
n = length(d);
for k = 2:n
  lambda = c(k-1)/d(k-1);
  d(k) = d(k) - lambda*e(k-1);
  c(k-1) = lambda;
end