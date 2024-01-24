function root = polyroots(a, tol)

if nargin == 1
  tol = 1.0e-6;
end

n = length(a) - 1;
root = zeros(n,1);
for i = 1:n
  x = my_laguerre(a, tol);
  if abs(imag(x)) < tol
    x = real(x);
  end
  root(i) = x;
  a = deflpoly(a,x);
end
