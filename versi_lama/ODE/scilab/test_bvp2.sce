exec("ode_RK4.sce",-1)
exec("bvp2_shoot_secant.sce",-1)

function dx = df661(t,x)
  dx(1) = x(2)
  dx(2) = (2*x(1) + 4*t*x(2))*x(1)
  // don't forget to transpose
  dx = dx'
endfunction

t0 = 0
tf = 1
x0 = 1/4
xf = 1/3 // initial/final times and positions

N = 100
tol = 1e-8
kmax = 10

[t,x] = bvp2_shoot_secant(df661,t0,tf,x0,xf,N,tol,kmax)

// analytic solution
xo = 1.0 ./ (4 - t.*t)

err = norm(x(:,1) - xo)/(N + 1)

//clf()
//plot(t,x(:,1),'b', t,xo,'r')
