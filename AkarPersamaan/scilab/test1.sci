clear

exec "bisection.sci"
exec "regula_falsi.sci"
exec "ridder.sci"
exec "newton_raphson.sci"
exec "secant.sci"

function f = fun1(x)
  f = x^3 - 10*x^2 + 5
endfunction

function df = dfun1(x)
  df = 3*x^2 - 20*x
endfunction

//root = bisection( fun1, 0.6, 0.8 )
//root = regula_falsi( fun1, 0.6, 0.8 )
//root = ridder( fun1, 0.6, 0.8 )
//root = newton_raphson( fun1, dfun1, 7.0, 7.1 )
//root = newton_raphson( fun1, dfun1, 0.5 )
root = secant( fun1, 0.5 )

function f = fun2(x)
  f = cos(x) - x^3
endfunction

//root = bisection( fun2, 0.0, 1.0 )
//root = regula_falsi( fun2, 0.0, 1.0 )
//root = ridder( fun2, 0.0, 1.0 )


function f = fun3(x)
  f = 1/((x-0.3)^3 + 0.01) - 1/( (x-0.8)^2 + 0.04 )
endfunction

//root = bisection( fun3, 0.5, 0.7 )
//root = regula_falsi( fun3, 0.5, 0.7 )
//root = ridder( fun3, 0.5, 0.7 )

exit
