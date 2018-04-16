exec("wave_1d.sce",-1)
exec("to_string.sce",-1)

function y = it0(x)
  omega = 2*%pi/5.0
  y = sin(omega*x) // + 0.25*sin(2*omega*x) + 0.25*sin(4*omega*x)
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
tf = 5.0
Nt = 100

dx = xf/Nx
dt = tf/Nt

[u,x,t] = wave_1d(a,xf,tf,it0,i1t0,bx0t,bxft,Nx,Nt)

for n = 1:Nt
  clf()
  plot(x,u(:,n))
  set(gca(),"data_bounds", [0 xf -1  1])
  xs2png(gcf(), "TEMP_wave_" + to_string(n) + ".png")
end