exec("heat_2d_ADI.sce",-1)

function z = it0(x,y)
  z = 0
endfunction

function z = bxyt(x,y,t)
  z = exp(y)*cos(x) - exp(x)*cos(y)
endfunction

a = 1e-4
D = [0 4 0 4]
T = 5000
Mx = 40
My = 40
N = 50

[u,x,y,t] = heat_2d_ADI(a,D,T,it0,bxyt,Mx,My,N)

mesh(x,y,u)
