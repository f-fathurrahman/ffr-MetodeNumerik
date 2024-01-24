exec("newton_raphson.sce", 0)

function f = fun1(x)
  f = x^3 - 10*x^2 + 5
endfunction

// derivative of fun1, required for newton_rapshon
function df = dfun1(x)
  df = 3*x^2 - 20*x
endfunction

root = newton_raphson( fun1, dfun1, 0.7 )

if getscilabmode() ~= "STD"
  quit()
end