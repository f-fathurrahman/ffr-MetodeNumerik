exec("to_string.sce",-1)
exec("heat_1d_euler_exp.sce",-1)
exec("heat_1d_euler_imp.sce",-1)
exec("heat_1d_CN.sce",-1)

// initial condition (function of x)
function T = it0(x)
  T = sin(%pi*x)
endfunction

// boundary condition (function of t)
function T = bx0(t)
  T = 0.0
endfunction

function T = bxf(t)
  T = 0.0
endfunction

function T = analytic_solution(x,t)
  T = sin(%pi*x)*exp(-%pi*%pi*t)
endfunction

function plot_to_png(u,x,t,prefix)
  Nt = length(t)-1
  for it = 1:Nt+1
    clf()
    plot(x,u(:,it))
    set(gca(), "data_bounds", [0,1,0,1])
    strt = "t = " + string(t(it))
    xstring(0.8,0.9,strt)
    xs2png(gcf(), prefix + to_string(it) + ".png")
    printf("Done output solution for t = %f\n", t(it))
  end
endfunction


a = 1

xf = 1
Nx = 25

T  = 0.1
Nt = 100

// Explicit Euler
[u1,x,t] = heat_1d_euler_exp( a, xf, T, it0, bx0, bxf, Nx, Nt )
plot_to_png(u1,x,t,"TEMP_exp_")

// Using implicit Euler method
[u2,x,t] = heat_1d_euler_imp( a, xf, T, it0, bx0, bxf, Nx, Nt )
plot_to_png(u2,x,t,"TEMP_imp_")

// Using Crank-Nicholson method
[u3,x,t] = heat_1d_CN( a, xf, T, it0, bx0, bxf, Nx, Nt )
plot_to_png(u3,x,t,"TEMP_CN_")


NxNt = Nx*Nt
u_analytic = analytic_solution(x,t)

//How far from the analytical solution?
err1 = norm((u1-u_analytic))/NxNt
err2 = norm((u2-u_analytic))/NxNt
err3 = norm((u3-u_analytic))/NxNt

printf("err1 = %f\n", err1)
printf("err2 = %f\n", err2)
printf("err3 = %f\n", err3)

if getscilabmode() ~= "STD"
  quit()
end
