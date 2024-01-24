// gaya: kN
// disp: cm
// This function is quite difficult to minimize because it has small curvature
function y = soal_02(x)
  w0 = 2.5
  E = 50000
  I = 30000
  L = 450
  C1 = -x^5 + 2*L^2*x^3 - L^4*x
  y = w0/(120*E*I*L)*C1
endfunction

function y = d_soal_02_x(x)
  w0 = 2.5
  E = 50000
  I = 30000
  L = 450  
  C1 = -5*x^4 + 6*L^2*x^2 - L^4
  y = w0/(120*E*I*L)*C1
endfunction

// without constant w0/(120*E*I*L)
function y = d_soal_02_x_v2(x)
  w0 = 2.5
  E = 50000
  I = 30000
  L = 450  
  C = w0/(120*E*I*L)
  printf("Constant = %20.18f\n", C)
  y = -5*x^4 + 6*L^2*x^2 - L^4
endfunction

function f = xfx_d_soal_02_x(x)
  w0 = 2.5
  E = 50000
  I = 30000
  L = 450  
  // C1 = -5*x^4 + 6*L^2*x^2 - L^4
  C = w0/(120*E*I*L)
  f = sqrt( abs((-5*x^4 - L^4)/(6*L^2)) ) // need to use absolute value
endfunction

function do_plot()
  x1 = 100.0
  x2 = 250.0
  Npoints = 100
  x = linspace(x1, x2, Npoints)
  y = zeros(1,Npoints)
  dydx = zeros(1,Npoints)
  for i = 1:Npoints
    y(i) = soal_02(x(i))
    dydx(i) = d_soal_02_x(x(i))
    printf("%18.10f %18.10f %18.10f\n", x(i), y(i), dydx(i))
  end
  //
  clf()
  plot(x, y)
  xgrid()
  xs2pdf( gcf(), "soal_02_y.pdf" )
  //
  clf()
  plot(x, dydx)
  xgrid()
  xs2pdf( gcf(), "soal_02_dydx.pdf" )  
endfunction

// do_plot()

exec("bisection.sce", -1)
root = bisection( d_soal_02_x, 190, 210 )
printf("At root bisection, derivative is %18.10f\n", d_soal_02_x(root))
printf("At root bisection, minimum deflection is %18.10f\n", soal_02(root))

exec("regula_falsi.sce", -1)
root = regula_falsi( d_soal_02_x, 190, 210 )
printf("At root regula_falsi, derivative is %18.10f\n", d_soal_02_x(root))
printf("At root regula_falsi, minimum deflection is %18.10f\n", soal_02(root))

exec("ridder.sce", -1)
root = ridder( d_soal_02_x, 190, 210 )
printf("At root ridder, derivative is %18.10f\n", d_soal_02_x(root))
printf("At root ridder, minimum deflection is %18.10f\n", soal_02(root))

exec("fixed_point.sce", -1)
root = fixed_point( xfx_d_soal_02_x, 190 )
printf("At root fixed_point, derivative is %18.10f\n", d_soal_02_x(root))
printf("At root fixed_point, minimum deflection is %18.10f\n", soal_02(root))

exec("secant.sce", -1)
root = secant( d_soal_02_x, 190 )
printf("At root secant, derivative is %18.10f\n", d_soal_02_x(root))
printf("At root secant, minimum deflection is %18.10f\n", soal_02(root))


if getscilabmode() ~= "STD"
  quit()
end
