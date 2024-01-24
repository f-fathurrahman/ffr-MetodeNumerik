function [A,perm] = LUdecPiv(A)
% LU decomposition of matrix A; returns A = [L\U]
% and the row permutation vector 'perm'.
% USAGE: [A,perm] = LUdecPiv(A)
n = size(A,1);
s = zeros(n,1);
perm = (1:n)';

% Set up scale factor array
for i = 1:n;
  s(i) = max(abs(A(i,1:n)));
end

% Exchange rows if necessary
for k = 1:n-1
  [Amax,p] = max(abs(A(k:n,k))./s(k:n));
  p = p + k - 1;
  if Amax < eps
    error('Matrix is singular')
  end
  if p ~= k
    s = swapRows(s,k,p);
    A = swapRows(A,k,p);
    perm = swapRows(perm,k,p);
  end
  % Elimination pass
  for i = k+1:n
    if A(i,k) ~= 0
      lambda = A(i,k)/A(k,k);
      A(i,k+1:n) = A(i,k+1:n) - lambda*A(k,k+1:n);
      A(i,k) = lambda;
    end
  end
end
