function b = deflpoly(a, r)
% Horner's deflation
n = length(a) - 1;
b = zeros(n,1);
b(1) = a(1);
for i = 2:n
  b(i) = a(i) + r*b(i-1);
end