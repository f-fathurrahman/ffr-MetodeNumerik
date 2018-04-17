exec("wave_1d.sce",-1)
exec("to_string.sce",-1)

function y = it0(x)
  x0 = 2.5
  dr2 = (x - x0)^2
  y = exp(-5*dr2)
endfunction

function y = i1t0(x)
  y = 0
endfunction

function y = bx0t(t)
  y = 0
endfunction

function y = bxft(t)
  y = 0
endfunction

a  = 1

xf = 5.0

Nx = 100
tf = 20.0
Nt = 500

dx = xf/Nx
dt = tf/Nt

[u,x,t] = wave_1d(a,xf,tf,it0,i1t0,bx0t,bxft,Nx,Nt)

for n = 1:Nt
  clf()
  plot(x,u(:,n))
  set(gca(),"data_bounds", [0 xf -1  1])
  xs2png(gcf(), "TEMP_wave_" + to_string(n) + ".png")
end