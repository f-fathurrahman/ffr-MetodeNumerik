exec("to_string.sce",-1)
exec("wave_2d.sce",-1)

function z = it0(x,y)
  kx = 2*%pi
  ky = 2*%pi
  z = sin(3*kx*x)*sin(ky*y)
endfunction

function z = i1t0(x,y)
  z = 0.0
endfunction

function z = bxyt(x,y,t)
  z = 0.0
endfunction

a = 0.25
D = [0 1 0 1]
T = 2
Mx = 40
My = 40
N = 100

[u,x,y,t] = wave_2d(a,D,T,it0,i1t0,bxyt,Mx,My,N)
