// Incremental search for a root of f(x)
// func: handle of function
// a, b: limits of search
// dx: search increment
// x1, x2: bounds on the smallest root in (a,b)
//         set to NaN if no root was detected
function [x1, x2] = rootsearch( func, a, b, dx )

  x1 = a
  f1 = func(x1)

  x2 = a + dx
  f2 = func(x2)

  while f1*f2 > 0.0
    if x1 >= b
      x1 = %nan
      x2 = %nan
      return
    end
    x1 = x2
    f1 = f2
    x2 = x1 + dx
    f2 = func(x2)
  end

endfunction
