A = [6 -4 1;
     -4 6 -4;
     1 -4 6]
b = [-14 22;
     36 -18;
      6   7]

[x1,det1] = gauss(A,b(:,1))

[x2,det2] = gauss(A,b(:,2))
