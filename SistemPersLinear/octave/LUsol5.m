function x = LUsol5(d,e,f,b)
% Solves A*x = b where A = [f\e\d\e\f] is the LU
% decomposition of the original pentadiagonal A.
% USAGE: x = LUsol5(d,e,f,b)
n = length(d);
b(2) = b(2) - e(1)*b(1);
% Forward substitution
for k = 3:n
  b(k) = b(k) - e(k-1)*b(k-1) - f(k-2)*b(k-2);
end
b(n) = b(n)/d(n);
% Back substitution
b(n-1) = b(n-1)/d(n-1) - e(n-1)*b(n);
for k = n-2:-1:1
  b(k) = b(k)/d(k) - e(k)*b(k+1) - f(k)*b(k+2);
end
x = b;
