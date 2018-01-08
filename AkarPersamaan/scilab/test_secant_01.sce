exec("secant.sce", 0)

function f = fun1(x)
  f = x^3 - 10*x^2 + 5
endfunction

root = secant( fun1, 0.7 )

if getscilabmode() ~= "STD"
  quit()
end