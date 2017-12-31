function x = LUsol3(c,d,e,b)
% Solves A*x = b where A = [c\d\e] is the LU
% decomposition of the original tridiagonal A.
% USAGE: x = LUsol3(c,d,e,b)
n = length(d);
for k = 2:n
  % Forward substitution
  b(k) = b(k) - c(k-1)*b(k-1);
end

b(n) = b(n)/d(n);

% Back substitution
for k = n-1:-1:1
  b(k) = (b(k) -e(k)*b(k+1))/d(k);
end
x = b;