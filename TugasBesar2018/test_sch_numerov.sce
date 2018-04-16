exec("sch_numerov.sce",-1)
exec("sch_RK4.sce",-1)

function V = HarmonicPot(x)
  V = 0.5 * x^2
endfunction

E = 0.5
xspan = [0 5.0]
N = 2000
h = (xspan(2) - xspan(1))/N
printf("h = %18.10f\n", h)

// for odd solution
//y0  = 0.0
//dy0 = 1.0

// for even solution
y0  = 1.0
dy0 = 0.0

[x,y,idx_div] = sch_numerov( E, HarmonicPot, xspan, y0 ,dy0 ,N )
plot(x,y,"r")
[x,y,idx_div] = sch_numerov( E-0.001, HarmonicPot, xspan, y0 ,dy0 ,N )
plot(x,y,"g")
[x,y,idx_div] = sch_numerov( E+0.001, HarmonicPot, xspan, y0 ,dy0 ,N )
plot(x,y,"b")
[x,y,idx_div] = sch_numerov( E-0.002, HarmonicPot, xspan, y0 ,dy0 ,N )
plot(x,y,"g--")
[x,y,idx_div] = sch_numerov( E+0.002, HarmonicPot, xspan, y0 ,dy0 ,N )
plot(x,y,"g--")
set( gca(),"data_bounds",[0.0 5.0 -1.0 1.0] )

// y(1) -> psi
// y(2) -> psi' -> f(1)
// f(2) -> y'(2) -> psi'' = 2*(V(x) - E)*y(1)
//function f = dy(E,x,y)
//  f(1) = y(2)
//  f(2) = 2*(HarmonicPot(x) - E)*y(1)
//  f = f'
//endfunction
//[x1,y1] = sch_RK4( E, dy, xspan, [y0 dy0], N )
//plot(x1,y1(:,1))
