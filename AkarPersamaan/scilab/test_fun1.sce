function [x, y] = myfct(a, b)
  x = a + b
  y = a - b
endfunction

[x, y] = myfct(3,2)

printf("x = %f, y = %f\n", x, y)

quit()
