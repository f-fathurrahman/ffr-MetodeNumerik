function V = Potential(x)
  if abs(x) > 1.0
    V = 1
  else
    V = 0
  end
endfunction

omega = 1
N = 200
L = 10.0  // domain is defined on [-L/2,L/2]

h = L/(N-1)

nabla2 = zeros(N,N)  // Laplacian
V = zeros(N)

x = -L/2 + [0:N-1]'*h

for i = 1:N-1
  nabla2(i,i) = -2
  nabla2(i+1,i) = 1
  nabla2(i,i+1) = 1
end
nabla2(N,N) = -2
nabla2(:,:) = nabla2/h/h

K = -0.5*nabla2
for i = 1:N
  V(i) = Potential(x(i))
end
H = K + diag(V)

[evecs,D] = spec(H)
evals = diag(D)
