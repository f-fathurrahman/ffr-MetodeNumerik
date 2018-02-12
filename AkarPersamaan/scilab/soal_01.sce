function f = soal_01(omega)
  Z = 75
  R = 225
  C = 0.6e-6
  L = 0.5
  if abs(omega) < %eps
    printf("omega is small: %18.10f\n", omega)
    return %nan
  end
  A = 1/R^2 + (omega*C - 1/(omega*L))^2
  f = 1/Z - sqrt(A)
endfunction

function f = d_soal_01(omega)
  f = omega
endfunction

function do_plot()
  Npoints = 100
  omega_1 = 50
  omega_2 = 200
  omega = linspace(omega_1, omega_2, Npoints)
  f = zeros(1,Npoints)
  for i = 1:Npoints
    f(i) = soal_01(omega(i))
    printf("%18.10f %18.10f\n", omega(i), f(i))
  end
  clf()
  plot(omega, f)
  xgrid()
  xs2pdf( gcf(), "soal_01.pdf" )
endfunction

// do_plot()

exec("bisection.sce", -1)
root = bisection( soal_01, 150, 170 )

exec("regula_falsi.sce", -1)
root = regula_falsi( soal_01, 150, 170 )

exec("ridder.sce", -1)
root = ridder( soal_01, 150, 170 )

exec("secant.sce", -1)
root = secant( soal_01, 160 )

if getscilabmode() ~= "STD"
  quit()
end
