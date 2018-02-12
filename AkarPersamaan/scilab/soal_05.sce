function f = soal_05(t)
  g = 9.81 // m*s^-2
  u = 2200 // m/s
  m0 = 160000 // kg
  q = 2680 // kg/s
  v = 1000 // m/s
  C = m0/(m0 - q*t)
  f = u*log(C) - g*t - v
endfunction

function v = eval_v(t)
  g = 9.81 // m*s^-2
  u = 2200 // m/s
  m0 = 160000 // kg
  q = 2680 // kg/s
  C = m0/(m0 - q*t)
  v = u*log(C) - g*t
endfunction

function do_plot()
  t1 = 1.0
  t2 = 40.0
  Npoints = 100
  t = linspace(t1, t2, Npoints)
  f = zeros(1,Npoints)
  for i = 1:Npoints
    f(i) = soal_05(t(i))
    printf("%18.10f %18.10f\n", t(i), f(i))
  end
  clf()
  plot(t, f)
  xgrid()
  xs2pdf( gcf(), "soal_05.pdf" )  
endfunction

//do_plot()

exec("bisection.sce", -1)
root = bisection( soal_05, 25, 26 )
printf("At root = %18.10f\n", eval_v(root))

exec("regula_falsi.sce", -1)
root = regula_falsi( soal_05, 25, 26 )
printf("At root = %18.10f\n", eval_v(root))

exec("ridder.sce", -1)
root = ridder( soal_05, 25, 26 )
printf("At root = %18.10f\n", eval_v(root))

if getscilabmode() ~= "STD"
  quit()
end
