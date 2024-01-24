function x = fex2_17(x,omega)
n = length(x);
x(1) = omega*(x(2) - x(n))/2 + (1-omega)*x(1);
for i = 2:n-1
  x(i) = omega*(x(i-1) + x(i+1))/2 + (1-omega)*x(i);
end
x(n) = omega *(1 - x(1) + x(n-1))/2 + (1-omega)*x(n);
