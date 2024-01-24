A = [6 -4 1;
     -4 6 -4;
     1 -4 6]
b = [-14 22;
     36 -18;
      6   7]

A_chol = choleskiDec(A)

x1 = choleskiSol(A_chol,b(:,1))
x2 = choleskiSol(A_chol,b(:,2))

disp('Test')
disp(A*x1 - b(:,1))
disp(A*x2 - b(:,2))


