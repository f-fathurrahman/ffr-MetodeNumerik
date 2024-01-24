function [t,x] = bvp2_shoot_secant(f,t0,tf,x0,xf,N,tol,kmax)
// To solve BVP2: [x1,x2]’ = f(t,x1,x2) with x1(t0) = x0, x1(tf) = xf
  
  if ~exists("kmax","local")
    kmax = 10
  end

  if ~exists("tol","local")
    tol = 1e-8
  end
  
  if ~exists("N","local")
    N = 100
  end

  dx0(1) = (xf - x0)/(tf - t0) // the initial guess of x'(t0)

  printf("Initial guess of dx(t0) = %f\n", dx0(1))

  [t,x] = ode_RK4(f,[t0 tf],[x0 dx0(1)],N)
  
//  clf()
//  plot(t,x(:,1))
//  //set( gca(), "data_bounds", [0,1,0.2,0.45] )
//  orig_data_bounds = get( gca(), "data_bounds" )
//  //disp(orig_data_bounds)
//  xs2png(gcf(),"TEMP_shoot_0.png")

  e(1) = x(N+1,1) - xf
  dx0(2) = dx0(1) - 0.1*sign(e(1))
  
  for k = 2:kmax-1

    [t,x] = ode_RK4(f,[t0 tf],[x0 dx0(k)],N)
    
    //clf()
    //plot(t,x(:,1))
    //set(gca(), "data_bounds", orig_data_bounds)
    //xs2png(gcf(),"TEMP_shoot_" + string(k) + '.png')
    
    //difference between the resulting final value and the target one
    e(k) = x(N+1,1) - xf

    printf("iter = %d, error = %18.10f\n", k, e(k))
  
    ddx = dx0(k) - dx0(k - 1)
  
    if abs(e(k))< tol | abs(ddx) < tol
      break
    end
  
    deddx = (e(k) - e(k - 1))/ddx // the gradient of mismatching error
  
    dx0(k+1) = dx0(k) - e(k)/deddx // move by secant method

  end

  Nk = length(dx0)
  printf("Final dx0 = %18.10f\n",dx0(Nk))

endfunction