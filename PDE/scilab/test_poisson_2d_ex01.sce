exec("poisson_2d.sce", -1)

function z = f(x,y)
  z = 0.0
endfunction

function z = g(x,y)
  z = 0.0
endfunction

function z = bx0(y)
  z = exp(y) - cos(y)
endfunction

function z = bxf(y)
  z = exp(y)*cos(4) - exp(4)*cos(y)
endfunction

function z = by0(x)
  z = cos(x) - exp(x)
endfunction

function z = byf(x)
  z = exp(4)*cos(x) - exp(x)*cos(4)
endfunction

x0 = 0
xf = 4
Mx = 20

y0 = 0
yf = 4
My = 20

D = [x0 xf y0 yf]
MaxIter = 500
tol = 1e-4

[U,x,y] = solve_poisson(f,g,bx0,bxf,by0,byf,D,Mx,My,tol,MaxIter)

plot3d(x,y,U)

if getscilabmode() ~= "STD"
  quit()
end
