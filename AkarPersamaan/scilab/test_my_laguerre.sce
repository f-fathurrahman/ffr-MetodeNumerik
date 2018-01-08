exec("evalpoly.sce", 0)
exec("deflpoly.sce", 0)
exec("my_laguerre.sce", 0)

a = [1, 3, 4, 5, 6, 7, 8]

my_laguerre(a, 1e-9)

if getscilabmode() ~= "STD"
  quit()
end