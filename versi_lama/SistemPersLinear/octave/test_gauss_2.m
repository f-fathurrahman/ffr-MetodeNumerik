A = vander(1:0.2:2);
b = [0 1 0 1 0 1]';

x_exact = [1250/3 -3125 9250 -13500 29128/3 -2751]';
[x,det1] = gauss(A,b)

norm(x-x_exact)

A*x

A*x_exact
