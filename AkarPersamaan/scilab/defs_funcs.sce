function f = fun1(x)
  f = x^3 - 10*x^2 + 5
endfunction

// derivative of fun1, required for newton_rapshon
function df = dfun1(x)
  df = 3*x^2 - 20*x
endfunction