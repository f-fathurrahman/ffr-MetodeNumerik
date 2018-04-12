exec("sch_numerov.sce",-1)
exec("sch_RK4.sce",-1)

function V = HarmonicPot(x)
  V = 0.5 * x^2
endfunction

// y(1) -> psi
// y(2) -> psi' -> f(1)
// f(2) -> y'(2) -> psi'' = 2*(V(x) - E)*y(1)
function f = dy(E,x,y)
  f(1) = y(2)
  f(2) = 2*(HarmonicPot(x) - E)*y(1)
  f = f'
endfunction

E = 3.5
xspan = [0 5]
N = 100

// for even solution
// y0  = 1.0
// dy0 = 0.0

// for odd solution
y0  = 0.0
dy0 = 1.0

[x,y,idx_div] = sch_numerov( E, HarmonicPot, xspan, y0 ,dy0 ,N )

[x1,y1] = sch_RK4( E, dy, xspan, [y0 dy0], N )
