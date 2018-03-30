exec("heat_1d_euler_exp.sce",-1)

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
  T = sin(pi*x)*exp(-%pi*%pi*t)
endfunction

a = 1 // the parameter of (E9.2-1)

xf =   1
Nt =  50
T  =   0.2
Nx = 500

[u1,x,t] = heat_1d_euler_exp( a, xf, T, it0, bx0, bxf, Nt, Nx )

for it = 1:Nt+1
  clf()
  plot(x,u1(:,it))
  set(gca(), "data_bounds", [0,1,0,1])
  xs2png( gcf(), "TEMP_" + string(it) + ".png" )
  printf("Done output solution for time %d\n", it)
end

//figure(1), clf, mesh(t,x,u1)
// [u2,x,t]=heat_imp(a,xf,T,it0,bx0,bxf,M,N); %converge unconditionally
// figure(2), clf, mesh(t,x,u2)
// [u3,x,t]=heat_CN(a,xf,T,it0,bx0,bxf,M,N); %converge unconditionally
// figure(3), clf, mesh(t,x,u3)
// MN=M*N;
// Uo= uo(x,t); aUo=abs(Uo)+eps; // values of true analytical solution
// %How far from the analytical solution?
// err1= norm((u1-Uo)./aUo)/MN
// err2= norm((u2-Uo)./aUo)/MN
// err3= norm((u3-Uo)./aUo)/MN

if getscilabmode() ~= "STD"
  quit()
end
