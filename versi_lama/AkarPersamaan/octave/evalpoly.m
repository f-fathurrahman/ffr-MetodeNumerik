function [p, dp, ddp] = evalpoly(a,x)
% evaluate polynomial
% p = a(1)*x^n + a(2)*x^(n-1) + ... + a(n+1)
% and its first two derivatives dp and ddp

n = length(a) - 1;

p = a(1);
dp = 0.0;
ddp = 0.0;

for i = 1:n
  ddp = ddp*x + 2.0*dp;
  dp = dp*x + p;
  p = p*x + a(i+1);
end
