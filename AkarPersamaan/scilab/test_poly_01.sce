exec("evalpoly.sce", 0)
exec("deflpoly.sce", 0)

a = [1, 3, 4, 5, 6, 7, 8]

[p, dp, d_dp] = evalpoly(a,1)

x1 = deflpoly(a, 1.0)

if getscilabmode() ~= "STD"
  quit()
end
