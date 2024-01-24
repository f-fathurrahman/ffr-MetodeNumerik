function f = soal_06(z)
  epsilon0 = 0.885e-12
  F = 0.3
  Q = 9.4e-6
  q = 2.4e-5
  R = 0.1
  c = 1 - z/sqrt(z^2 + R^2)
  f = Q*q*c/(2*epsilon0) - F
endfunction

function F = eval_F(z)
  epsilon0 = 0.885e-12
  Q = 9.4e-6
  q = 2.4e-5
  R = 0.1
  c = 1 - z/sqrt(z^2 + R^2)
  F = Q*q*c/(2*epsilon0)
endfunction

function do_plot()
  z1 = 0.1
  z2 = 2.0
  Npoints = 100
  z = linspace(z1, z2, Npoints)
  f = zeros(1,Npoints)
  for i = 1:Npoints
    f(i) = soal_06(z(i))
    printf("%18.10f %18.10f\n", z(i), f(i))
  end
  clf()
  plot(z, f)
  xgrid()
  xs2pdf( gcf(), "soal_06.pdf" )  
endfunction

// do_plot()

exec("bisection.sce", -1)
root = bisection( soal_06, 1.4, 1.6 )
printf("At root = %18.10f\n", eval_F(root))

exec("regula_falsi.sce", -1)
root = regula_falsi( soal_06, 1.4, 1.6 )
printf("At root = %18.10f\n", eval_F(root))

exec("ridder.sce", -1)
root = ridder( soal_06, 1.4, 1.6 )
printf("At root = %18.10f\n", eval_F(root))

if getscilabmode() ~= "STD"
  quit()
end
