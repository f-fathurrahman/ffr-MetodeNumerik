exec("to_string.sce",-1)
exec("wave_2d.sce",-1)

function z = it0(x,y)
  z = 0.1*sin(%pi*x)*sin(%pi*y/2)
endfunction

function z = i1t0(x,y)
  z = 0.0
endfunction

function z = bxyt(x,y,t)
  z = 0.0
endfunction

a = 0.25
D = [0 2 0 2]
T = 2
Mx = 40
My = 40
N = 40

[u,x,y,t] = wave_2d(a,D,T,it0,i1t0,bxyt,Mx,My,N)
