function f = soal_10(omega)
  m = 2500 // kg
  k = 300000 // N/m
  c = 36e3 // Ns/m
  XperY = 0.4

  num = omega^2*c + k^2
  denum = (k - m*omega^2) + (omega*c)^2
  f = sqrt(num/denum) - XperY
endfunction

function do_plot()
  omega1 = 10
  omega2 = 30
  Npoints = 100
  omega = linspace(omega1, omega2, Npoints)
  f = zeros(1,Npoints)
  for i = 1:Npoints
    f(i) = soal_10(omega(i))
    printf("%18.10f %18.10f\n", omega(i), f(i))
  end  
endfunction

//do_plot()

exec("bisection.sce", -1)
root = bisection( soal_10, 20.7070707071, 20.9090909091 )

exec("regula_falsi.sce", -1)
root = regula_falsi( soal_10, 20.7070707071, 20.9090909091 )

exec("ridder.sce", -1)
root = ridder( soal_10, 20.7070707071, 20.9090909091 )

if getscilabmode() ~= "STD"
  quit()
end
