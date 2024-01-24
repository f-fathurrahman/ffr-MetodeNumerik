function [u,x,y,t] = wave_2d(a,D,T,it0,i1t0,bxyt,Mx,My,N)
 // solve a(u_xx +u_yy)=u_tt for D(1)<=x<=D(2),D(3)<=y<=D(4), 0<=t<=T
 //  Initial Condition: u(x,y,0) = it0(x,y), u_t(x,y,0) = i1t0(x,y)
 //  Boundary Condition: u(x,y,t) = bxyt(x,y,t) for (x,y)cBoundary
 //  Mx/My = # of subintervals along x/y-axis
 //  N = # of subintervals along t-axis
  dx = ( D(2) - D(1) )/Mx
  x  = D(1) + [0:Mx]*dx

  dy = ( D(4) - D(3) )/My
  y  = D(3) + [0:My]'*dy
  dt = T/N
  t  = [0:N]*dt
  
  // Initialization
  u  = zeros(My+1,Mx+1)
  ut = zeros(My+1,Mx+1)
  
  for j=2:Mx
    for i=2:My
      u(i,j)  = it0(x(j),y(i))
      ut(i,j) = i1t0(x(j),y(i))
    end
  end

  adt2 = a*dt*dt
  rx = adt2/(dx*dx)
  ry = adt2/(dy*dy)
  rxy1 = 1 - rx - ry
  rxy2 = rxy1*2
  
  u_1 = u
  
  for k = 0:N
    
    t = k*dt
    
    // should not have spatial (x and y) dependence
    for i = 1:My+1
      u(i,1) = bxyt(x(1),y(i),t)
      u(i,Mx+1) = bxyt(x(Mx+1),y(i),t)
    end

    for j = 1:Mx+1
      u(1,j) = bxyt(x(j),y(1),t)
      u(My+1,j)= bxyt(x(j),y(My+1),t)
    end
    
    if k==0 // starting condition
  
      for i=2:My
        for j=2:Mx
          u(i,j) = 0.5*(rx*(u_1(i,j-1) + u_1(i,j+1)) + ...
                   ry*(u_1(i-1,j) + u_1(i+1,j))) + rxy1*u(i,j) + dt*ut(i,j)
        end
      end

    else // propagation

      for i=2:My
        for j=2:Mx
          u(i,j) = rx*(u_1(i,j-1) + u_1(i,j+1)) + ...
                   ry*(u_1(i-1,j) + u_1(i+1,j)) + rxy2*u(i,j) - u_2(i,j)
        end
      end
    
    end
    
    u_2 = u_1
    u_1 = u
    
    clf()
    surf(x,y,u)
    set(gca(), "data_bounds", [0 1 0 1 -1 1])
    set(gcf(),"color_map",jetcolormap(32))
    f = gcf()
    f.figure_size = [800,1200]
    fc = f.children
    fc.rotation_angles = [77.75,-161.75]
    if k == 0
      xs2pdf(gcf(),"TEMP_wave2d_" + to_string(k) + ".pdf")
      //xs2png(gcf(),"TEMP_wave2d_" + to_string(k) + ".png")
    else
      //xs2png(gcf(),"TEMP_wave2d_" + to_string(k) + ".png")
      xs2pdf(gcf(),"TEMP_wave2d_" + to_string(k) + ".pdf")
    end

  end

endfunction

