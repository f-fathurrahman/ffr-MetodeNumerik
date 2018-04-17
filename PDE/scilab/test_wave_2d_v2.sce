exec("to_string.sce",-1)
exec("wave_2d.sce",-1)

function z = it0(x,y)
  cx = 0.5
  cy = 0.5
  dr2 = (x-cx)^2 + (y-cy)^2
  z = exp(-50.0*dr2)
endfunction

function z = i1t0(x,y)
  z = 0.0
endfunction

function z = bxyt(x,y,t)
  z = 0.0
endfunction

a = 0.25
D = [0 1 0 1]
T = 4
Mx = 40
My = 40
N = 200

[u,x,y,t] = wave_2d(a,D,T,it0,i1t0,bxyt,Mx,My,N)
