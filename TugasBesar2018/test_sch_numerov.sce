exec("sch_numerov.sce",-1)

function V = HarmonicPot(x)
  V = 0.5 * x^2
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

[x,y,div_idx] = sch_numerov(E,HarmonicPot,xspan,y0,dy0,N)