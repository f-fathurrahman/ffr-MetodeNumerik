exec("regula_falsi.sce", mode=0)

function f = fun1(x)
  f = x^3 - 10*x^2 + 5
endfunction

root = regula_falsi( fun1, 0.8, 0.6 )
printf("Final root = %18.10f\n", root)

quit()