exec("poisson2d.sce",-1)

function u = func(x,y)
  u = cos(x+y) - cos(x-y)
endfunction

function u = gauss(x,y)
  r2 = x**2 + y**2
  u = exp(-0.5*r2)
endfunction

Nx = 50
Ny = 50
x = linspace(-%pi,%pi,Nx)
y = linspace(-%pi,%pi,Ny)

u0 = zeros(Nx,Ny)
// set BC
u0(1,:) = 0.0
u0(:,1) = 0.0
u0(Nx,:) = 0.0
u0(:,Ny) = 0.0

// calculate array for RHS
f = zeros(Nx,Ny)
for j = 1:Ny
  for i = 1:Nx
    f(i,j) = func( x(i), y(j) )
  end
end

u = poisson2d(u0,x,y,Nx,Ny,1e-5,f)

surf(x,y,u)
set(gcf(),"color_map",jetcolormap(32))
colorbar(min(u),max(u))
xs2pdf(gcf(),"poisson2d.pdf")