exec("sch_numerov.sce",-1)
exec("sch_shoot.sce",-1)

function V = HarmonicPot(x)
  V = 0.5 * x^2
endfunction

// for odd solution
//y0  = 0.0
//dy0 = 1.0

// for even solution
y0  = 1.0
dy0 = 0.0

xspan = [0 10]
N = 2000
//idx_check = 60
EPS = 1e-6

dE_search = 0.1
Emin_search = 0.1
Emax_search = 15.0

E1 = Emin_search
E2 = Emin_search + dE_search

while E2 < Emax_search
  //printf("Searching root within [%f %f]\n", E1, E2)
  [E,is_exist] = sch_shoot(E1,E2,HarmonicPot,xspan,y0,dy0,N,EPS)
  E1 = E2
  E2 = E2 + dE_search
end