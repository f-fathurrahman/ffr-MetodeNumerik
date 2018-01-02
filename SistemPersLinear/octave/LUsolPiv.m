function x = LUsolPiv(A,b,perm)
% Solves L*U*b = x, where A contains row-wise
% permutation of L and U in the form A = [L\U].
% Vector ’perm’ holds the row permutation data.
% USAGE: x = LUsolPiv(A,b,perm)

% Rearrange b, store it in x
if size(b) > 1
  b = b'
end

n = size(A,1);
x = b;

for i = 1:n
  x(i) = b(perm(i));
end

% Forward and back substitution
for k = 2:n
  x(k) = x(k) - A(k,1:k-1)*x(1:k-1);
end

for k = n:-1:1
  x(k) = (x(k) - A(k,k+1:n)*x(k+1:n))/A(k,k);
end
