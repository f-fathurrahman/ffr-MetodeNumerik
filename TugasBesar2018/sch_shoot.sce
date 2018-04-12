function [E,is_exist] = sch_shoot( E1, E2, V, xspan, y0, dy0, N, EPS )
  MaxIter = 100

  [x,y,idx_div] = sch_numerov(E1, V, xspan, y0, dy0, N)
  f1 = y(idx_div)

  [x,y,idx_div] = sch_numerov(E2, V, xspan, y0, dy0, N)
  f2 = y(idx_div)

  // use bisection
  if f1*f2 < 0

    is_exist = %T
    
    iter = 1
    E_old = 1e10 // unlikely to be a solution

    while iter <= MaxIter+1

      E = 0.5*(E1 + E2)
      
      deltaE = abs(E_old - E)
      if deltaE < EPS
        printf("Root is found: E = %18.10f, deltaE = %18.10f\n", E, deltaE)
        break
      end
      
      [x,y,idx_div] = sch_numerov(E, V, xspan, y0, dy0, N)
      f = y(idx_div)

      if f1*f > 0
        // no root within [f1,f], change E1 to E
        E_old = E
        E1 = E
      else
        E_old = E
        E2 = E
      end

      iter = iter + 1

    end

    // print warning message if no root was found
    //if iter >= MaxIter
    //  printf("sch_shoot: No root was found within %f and %f\n", E1, E2)
    //end

  else

    //printf("Root is not found in interval %f and %f\n", E1, E2)
    E = 1e10
    is_exist = %F

  end


endfunction