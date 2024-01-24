exec("evalpoly.sce", 0)
exec("deflpoly.sce", 0)
exec("my_laguerre.sce", 0)
exec("polyroots.sce", 0)

a = [1, 3, 4, 5, 6, 7, 8];

r0 = polyroots(a)

r0_sci = roots(a)
n = length(r0)
printf("\nResult of built-in roots function of Scilab:\n")
for i = 1:n
  printf("Root %d = %18.10f %18.10f\n", i, real(r0_sci(i)), imag(r0_sci(i)))
end

if getscilabmode() ~= "STD"
  quit()
end
