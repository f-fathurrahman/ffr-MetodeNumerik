// gaya: kN
// disp: cm
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

exec("bisection.sce", -1)
root = bisection( d_soal_02_x, 190, 210 )

// clf()
// plot(x, y)
// xs2pdf( gcf(), "soal_02_y.pdf" )

// clf()
// plot(x, dydx)
// xs2pdf( gcf(), "soal_02_dydx.pdf" )

if getscilabmode() ~= "STD"
  quit()
end
