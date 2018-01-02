A = [6 -4 1;
     -4 6 -4;
     1 -4 6]
b = [-14 22;
     36 -18;
      6   7]

A_LU = LUdec(A)

x1 = LUsol(A_LU,b(:,1))

x2 = LUsol(A_LU,b(:,2))
