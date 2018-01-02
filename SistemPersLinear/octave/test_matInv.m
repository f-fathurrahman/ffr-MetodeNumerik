A = [ 0.6 -0.4 1.0;
     -0.3  0.2 0.5;
      0.6 -1.0 0.5];
Ainv = matInv(A)
check = A*Ainv