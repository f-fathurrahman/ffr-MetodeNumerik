function f = soal_01(Z, R, C, L, omega)
  A = 1/R^2 + (omega*C - 1/(omega*L))^2
  f = 1/Z - sqrt(A)
endfunction

function f = soal_01_ex1( omega )
  Z = 75
  R = 225
  C = 0.6e-6
  L = 0.5
  f = soal_01(Z, R, C, L, omega)
endfunction

Npoints = 100
omega_1 = 50
omega_2 = 200
omega = linspace(omega_1, omega_2, Npoints)
disp(size(omega))
N = disp(length(omega))
f = zeros(1,Npoints)
for i = 1:Npoints
  f(i) = soal_01_ex1(omega(i))
  printf("%18.10f %18.10f\n", omega(i), f(i))
end

clf()
plot(omega, f)
xs2pdf( gcf(), "soal_01.pdf" )

exec("bisection.sce", -1)
root = bisection( soal_01_ex1, 150, 170 )


if getscilabmode() ~= "STD"
  quit()
end
