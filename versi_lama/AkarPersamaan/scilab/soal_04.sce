 V = 3 // m^3

function f = soal_04(v)
  // Constants
  Ru = 0.518 // kJ/(kg K)
  // methane
  Pc = 4600 // kPa
  Tc = 191 // Kelvin
  //
  a = 0.427 * Ru^2 * Tc^2.5 / Pc
  b = 0.0866 * Ru * Tc / Pc
  T = -40 + 273 // K
  P = 65000 // kPa  
  //
  denum1 = v - b
  denum2 = v*(v + b)*sqrt(T)
  f = Ru*T*denum2 - a*denum1 - P*denum1*denum2
endfunction

function P = eval_P(v)
  // Constants
  Ru = 0.518 // kJ/(kg K)
  // methane
  Pc = 4600 // kPa
  Tc = 191 // Kelvin
  //
  a = 0.427 * Ru^2 * Tc^2.5 / Pc
  b = 0.0866 * Ru * Tc / Pc
  T = -40 + 273 // K
  P = 65000 // kPa  
  //  
  denum1 = v - b
  denum2 = v*(v + b)*sqrt(T)
  P = Ru*T/denum1 - a/denum2
endfunction

function do_plot()
  v1 = 0.001
  v2 = 0.005
  Npoints = 100
  v = linspace(v1, v2, Npoints)
  f = zeros(1,Npoints)
  for i = 1:Npoints
    f(i) = soal_04(v(i))
    printf("%18.10f %18.10f\n", v(i), f(i))
  end
  clf()
  plot(v, f)
  xgrid()
  xs2pdf( gcf(), "soal_04.pdf" )  
endfunction

// do_plot()

exec("bisection.sce", -1)
root = bisection( soal_04, 0.0025, 0.0030 )
printf("At root = %18.10f\n", eval_P(root))

exec("regula_falsi.sce", -1)
root = regula_falsi( soal_04, 0.0025, 0.0030 )
printf("At root = %18.10f\n", eval_P(root))

exec("ridder.sce", -1)
root = ridder( soal_04, 0.0025, 0.0030 )
printf("At root = %18.10f\n", eval_P(root))

if getscilabmode() ~= "STD"
  quit()
end
