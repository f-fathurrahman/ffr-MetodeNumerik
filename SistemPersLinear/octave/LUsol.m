function x = LUsol(A,b)
% Solves L*U*b = x, where A contains both L and U;
% that is, A has the form [L\U].
% USAGE: x = LUsol(A,b)
if size(b,2) > 1
  b = b';
end

n = length(b);

for k = 2:n
  b(k) = b(k) - A(k,1:k-1)*b(1:k-1);
end

for k = n:-1:1
  b(k) = (b(k) - A(k,k+1:n)*b(k+1:n))/A(k,k);
end

x = b;