exec("ode_euler.sce",-1)
exec("ode_heun.sce",-1)
exec("ode_RK4.sce",-1)

function x = df61(t,y)
  x = -y + 1
endfunction

a = 1
r = 1
tspan = [0 2]
y0 = 0
t = tspan(1) + [0:100]*(tspan(2)-tspan(1))/100
yt = 1 - exp(-a*t)

N = 4

[t1,ye] = ode_euler(df61,tspan,y0,N)
[t2,yh] = ode_heun(df61,tspan,y0,N)
[t3,yr] = ode_RK4(df61,tspan,y0,N)

clf()
plot(t,yt,"k")
plot(t1,ye,"b")
plot(t2,yh,"r")
plot(t3,yr,"g")

legend("exact","ode_euler","ode_heun", "ode_RK4")


