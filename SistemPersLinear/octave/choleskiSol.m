function x = choleskiSol(L,b)
% Solves [L][L’]{x} = {b}
% USAGE: x = choleskiSol(L,b)
n = length(b);

% {b} must be column vector
if size(b,2) > 1
  b = b'
end

% Solution of [L]{y} = {b}
for k = 1:n
  b(k) = (b(k) - dot(L(k,1:k-1),b(1:k-1)'))/L(k,k);
end

% Solution of {L}'{x} = {y}
for k = n:-1:1
  b(k) = (b(k) - dot(L(k+1:n,k),b(k+1:n)))/L(k,k);
end

x = b;