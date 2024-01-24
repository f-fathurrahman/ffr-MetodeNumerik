function Ainv = matInv(A)
% Inverts martix A with LU decomposition.
% USAGE: Ainv = matInv(A)

n = size(A,1);
Ainv = eye(n); % Store RHS vectors in Ainv.

[A,perm] = LUdecPiv(A); % Decompose A.

% Solve for each RHS vector and store results in Ainv
% replacing the corresponding RHS vector.
for i = 1:n
  Ainv(:,i) = LUsolPiv(A,Ainv(:,i),perm);
end
