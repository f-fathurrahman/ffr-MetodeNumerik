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
  T = sin(pi*x)*exp(-%pi*%pi*t)
endfunction

a = 1

xf = 1
Nx = 25

T  = 0.1
Nt = 100

// Explicit Euler

// [u1,x,t] = heat_1d_euler_exp( a, xf, T, it0, bx0, bxf, Nx, Nt )

//for it = 1:Nt+1
//  clf()
//  plot(x,u1(:,it))
//  set(gca(), "data_bounds", [0,1,0,1])
//  xs2png( gcf(), "TEMP_exp_" + to_string(it) + ".png" )
//  printf("Done output solution for t = %f\n", t(it))
//end

// Using implicit Euler method

// [u2,x,t] = heat_1d_euler_imp( a, xf, T, it0, bx0, bxf, Nx, Nt )

//for it = 1:Nt+1
//  clf()
//  plot(x,u2(:,it))
//  set(gca(), "data_bounds", [0,1,0,1])
//  xs2png( gcf(), "TEMP_imp_" + to_string(it) + ".png" )
//  printf("Done output solution for t = %f\n", t(it))
//end

// Using Crank-Nicholson method

[u3,x,t] = heat_1d_CN( a, xf, T, it0, bx0, bxf, Nx, Nt )

for it = 1:Nt+1
  clf()
  plot(x,u3(:,it))
  set(gca(), "data_bounds", [0,1,0,1])
  xs2png( gcf(), "TEMP_CN_" + to_string(it) + ".png" )
  printf("Done output solution for t = %f\n", t(it))
end

// MN=M*N;
// Uo= uo(x,t); aUo=abs(Uo)+eps; // values of true analytical solution
// %How far from the analytical solution?
// err1= norm((u1-Uo)./aUo)/MN
// err2= norm((u2-Uo)./aUo)/MN
// err3= norm((u3-Uo)./aUo)/MN

if getscilabmode() ~= "STD"
  quit()
end
