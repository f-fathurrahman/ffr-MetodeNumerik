function L = choleskiDec(A)
% Computes L in Choleski’s decomposition A = LL’.
% USAGE: L = choleski(A)
n = size(A,1);
for j = 1:n
  temp = A(j,j) - dot(A(j,1:j-1),A(j,1:j-1));
  if temp < 0.0
    error('Matrix is not positive definite')
  end
  A(j,j) = sqrt(temp);
  for i = j+1:n
    A(i,j) = (A(i,j) - dot(A(i,1:j-1),A(j,1:j-1)))/A(j,j);
  end
end
% Return a new matrix formed by extracting the lower ('tril') or
% upper ('triu') triangular part of the matrix A, and setting all
% other elements to zero.
L = tril(A);
