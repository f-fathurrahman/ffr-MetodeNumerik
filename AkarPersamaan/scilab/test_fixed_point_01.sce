function f = fun1(x)
  f = x^3 - 10*x^2 + 5
endfunction

function f = fun2(x)
  // original f(x) = x^3 - 10*x^2 + 5 = 0
  // modified to x = sqrt((x^3 + 5)/10)
  f = sqrt((x^3 + 5)/10)
endfunction

exec("fixed_point.sce", 0)
root = fixed_point( fun2, 0.7 )
printf("At root, f = %18.10f\n", fun1(root))

if getscilabmode() ~= "STD"
  quit()
end
