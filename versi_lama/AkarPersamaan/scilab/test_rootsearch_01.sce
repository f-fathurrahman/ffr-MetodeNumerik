function y = func1(x)
  y = x^3 - 10*x^2 + 5
endfunction

exec("rootsearch.sce")

dx = 0.1
a = 0
b = 1

[x1,x2] = rootsearch( func1, a, b, dx )
printf("x1 = %f, x2 = %f\n", x1, x2)
