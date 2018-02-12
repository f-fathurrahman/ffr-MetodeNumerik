function f = soal_08(omega)
  L = 11e-3
  C = 8e-6
  R = 1000
  denum = sqrt( (1 - omega^2*L*C)^2 + (omega*R*C)^2 )
  G = 0.87
  f = omega*R*C/denum - G
endfunction

function G = eval_G(omega)
  L = 11e-3
  C = 8e-6
  R = 1000
  denum = sqrt( (1 - omega^2*L*C)^2 + (omega*R*C)^2 )
  G = omega*R*C/denum
endfunction

function do_plot_v1()
  omega1 = 50
  omega2 = 300
  Npoints = 100
  omega = linspace(omega1, omega2, Npoints)
  f = zeros(1,Npoints)
  for i = 1:Npoints
    f(i) = soal_08(omega(i))
    printf("%18.10f %18.10f\n", omega(i), f(i))
  end
endfunction

function do_plot_v2()
  omega1 = 50000
  omega2 = 52000
  Npoints = 100
  omega = linspace(omega1, omega2, Npoints)
  f = zeros(1,Npoints)
  for i = 1:Npoints
    f(i) = soal_08(omega(i))
    printf("%18.10f %18.10f\n", omega(i), f(i))
  end
endfunction

//do_plot_v1()
//do_plot_v2()

exec("bisection.sce", -1)
root1 = bisection( soal_08, 219.1919191919, 221.7171717172 )
printf("At root1 = %18.10f\n", eval_G(root1))
//
root2 = bisection( soal_08, 51737.3737373737, 51757.5757575758 )
printf("At root2 = %18.10f\n", eval_G(root2))

exec("bisection.sce", -1)
root1 = bisection( soal_08, 219.1919191919, 221.7171717172 )
printf("At root1 = %18.10f\n", eval_G(root1))
//
root2 = bisection( soal_08, 51737.3737373737, 51757.5757575758 )
printf("At root2 = %18.10f\n", eval_G(root2))


exec("regula_falsi.sce", -1)
root1 = regula_falsi( soal_08, 219.1919191919, 221.7171717172 )
printf("At root1 = %18.10f\n", eval_G(root1))
//
root2 = regula_falsi( soal_08, 51737.3737373737, 51757.5757575758 )
printf("At root2 = %18.10f\n", eval_G(root2))

exec("ridder.sce", -1)
root1 = ridder( soal_08, 219.1919191919, 221.7171717172 )
printf("At root1 = %18.10f\n", eval_G(root1))
//
root2 = ridder( soal_08, 51737.3737373737, 51757.5757575758 )
printf("At root2 = %18.10f\n", eval_G(root2))

if getscilabmode() ~= "STD"
  quit()
end
