exec("ode_RK4.sce",-1)

// parameters and constants
global g
global L
global k

g = 9.81
L = 1.0
k = 0.0

// Initial conditions
// u0 will be varied later when calling ode_RK4
du0 = 0.0

t0 = 0.0
tmax = 20.0
delta_t = 0.01
N = (tmax - t0)/delta_t
printf("N = %d\n",N)

// k*u(2) is velocity dependent damping term
function f = soal_02(t,u)
  // make sure that we use the global variables
  global g
  global L
  global k
  //  
  f(1) =  u(2)
  f(2) = -g/L*sin(u(1)) - k*u(2)
  f = f'
endfunction

[t,u1] = ode_RK4( soal_02, [t0 tmax], [0.1 du0], N)
[t,u2] = ode_RK4( soal_02, [t0 tmax], [0.5*%pi du0], N)
[t,u3] = ode_RK4( soal_02, [t0 tmax], [3.1 du0], N)

//[t,u1] = ode_RK4( soal_02, [t0 tmax], [0.1 du0], N)
//[t,u2] = ode_RK4( soal_02, [t0 tmax], [0.1 du0], N)
//[t,u3] = ode_RK4( soal_02, [t0 tmax], [0.1 du0], N)

clf()
plot(t,u1(:,1),'r')
plot(t,u2(:,1),'g')
plot(t,u3(:,1),'b')
