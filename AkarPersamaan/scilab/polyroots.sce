function root = polyroots(a, tol)

  if ~exists("tol", "local")
    tol = 1e-6
  end

  n = length(a) - 1
  root = zeros(n,1)

  printf("\nFinding roots of polynomial:\n")
  for i=1:n
    printf("Power: %5d, coef = %18.10f\n", n - i, a(i))
  end
  printf("Tolerance = %e\n", tol)
  printf("\n")

  for i = 1:n
    x = my_laguerre(a, tol)
    is_complex_root = 1
    if abs(imag(x)) < tol
      x = real(x)
      is_complex_root = 0
    end
    root(i) = x
    //
    if is_complex_root
      printf("Found complex root %5d: %18.10f + %18.10fi\n", i, real(root(i)), imag(root(i)))
    else
      printf("Found real    root %5d: %18.10f\n", i, root(i))
    end
    //
    a = deflpoly(a,x)
  end

endfunction
