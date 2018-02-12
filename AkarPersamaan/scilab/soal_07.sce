function f = soal_07(h)
  r = 0.45 // 0.5*0.9
  mf = 70
  rho = 1030
  d = 2*r - h
  V_celup = (1.0/3.0) * %pi * d^2 * (3*r - d)
  f = rho*V_celup - mf
endfunction

function mf = eval_mf(h)
  r = 0.45 // 0.5*0.9
  mf = 70
  rho = 1030
  d = 2*r - h
  V_celup = (1.0/3.0) * %pi * d^2 * (3*r - d)
  mf = rho*V_celup
endfunction

function do_plot()
  h1 = 0.0
  h2 = 0.9
  Npoints = 100
  h = linspace(h1, h2, Npoints)
  f = zeros(1,Npoints)
  for i = 1:Npoints
    f(i) = soal_07(h(i))
    printf("%18.10f %18.10f\n", h(i), f(i))    
  end
  clf()
  plot(h, f)
  xgrid()
  xs2pdf( gcf(), "soal_07.pdf" )    
endfunction

// do_plot()

exec("bisection.sce", -1)
root = bisection( soal_07, 0.6545454545, 0.6636363636 )
printf("At root = %18.10f\n", eval_mf(root))

exec("regula_falsi.sce", -1)
root = regula_falsi( soal_07, 0.6545454545, 0.6636363636 )
printf("At root = %18.10f\n", eval_mf(root))

exec("ridder.sce", -1)
root = ridder( soal_07, 0.6545454545, 0.6636363636 )
printf("At ridder = %18.10f\n", eval_mf(root))

if getscilabmode() ~= "STD"
  quit()
end
