function [d,e,f] = LUdec5(d,e,f)
% LU decomposition of pentadiagonal matrix A = [f\e\d\e\f].
% USAGE: [d,e,f] = LUdec5(d,e,f)
n = length(d);

for k = 1:n-2
  lambda = e(k)/d(k);
  d(k+1) = d(k+1) - lambda*e(k);
  e(k+1) = e(k+1) - lambda*f(k);
  e(k) = lambda;
  lambda = f(k)/d(k);
  d(k+2) = d(k+2) - lambda*f(k);
  f(k) = lambda;
end
lambda = e(n-1)/d(n-1);
d(n) = d(n) - lambda*e(n-1);
e(n-1) = lambda;