function f = soal_03(x)
  K = 0.4
  pt = 3.5
  f = K - x/(1-x)*sqrt(2*pt/(2+x))
endfunction

function K = eval_K(x)
  pt = 3.5
  K = x/(1-x)*sqrt(2*pt/(2+x))
endfunction

x1 = 0.1
x2 = 0.3
Npoints = 100
x = linspace(x1,x2,Npoints)
f = zeros(1,Npoints)
for i = 1:100
  f(i) = soal_03(x(i))
  K = eval_K(x(i))
  printf("%18.10f %18.10f %18.10f\n", x(i), f(i), K)
end

// clf()
// plot(x, f)
// xs2pdf( gcf(), "soal_03.pdf" )

exec("bisection.sce", -1)
root = bisection( soal_03, x1, x2 )

K = eval_K(root)
printf("At root value of K is %18.10f\n", K)

if getscilabmode() ~= "STD"
  quit()
end

