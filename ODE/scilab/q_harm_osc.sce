function V = Potential(omega,x)
  //y = - 1/(abs(x));
  V = 0.5 * omega^2 * x.^2
endfunction

omega = 1
N = 200
L = 10.0  // domain is defined on [-L/2,L/2]

h = L/(N-1)

nabla2 = zeros(N,N)  // Laplacian
V = zeros(N)

t = -L/2 + [0:N-1]'*h

for i = 1:N-1
  nabla2(i,i) = -2
  nabla2(i+1,i) = 1
  nabla2(i,i+1) = 1
end
nabla2(N,N) = -2
nabla2(:,:) = nabla2/h/h

K = -0.5*nabla2
V = Potential(omega,t)
H = K + diag(V)

[evecs,D] = spec(H)
evals = diag(D)
