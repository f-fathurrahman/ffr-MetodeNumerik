A = [6 -4 1;
     -4 6 -4;
     1 -4 6]
b = [-14 22;
     36 -18;
      6   7]

[A_LU,perm] = LUdecPiv(A)

x1 = LUsolPiv(A_LU,b(:,1),perm)

x2 = LUsolPiv(A_LU,b(:,2),perm)
