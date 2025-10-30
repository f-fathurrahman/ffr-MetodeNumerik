! Start of Main Program
! Example 6.6
! ADI Method

    program main

    USE precisions
    USE DFPORT

! Declaration of variables
    IMPLICIT NONE
    INTEGER(int_p), PARAMETER :: nc = 300, mc = 60
    INTEGER(int_p), PARAMETER :: max_iter = 1000, max_times = 20
    INTEGER(int_p) :: i,j,iter,it
    REAL(real_p), PARAMETER :: one = 1.0D0,two = 2.0d0, &
                             & zero = 0.0d0,half = 0.5d0
    REAL(real_p), PARAMETER :: length = 0.5D0, height = 0.1d0, rho = 1000.0d0, &
                               JT0 = 10000.0d0, kap = 0.613d0, cp = 4200.0d0, &
                               phi_in = 20.0d0, dt = 10.0d0
    REAL(real_p) :: tol = 1.0d-4
    REAL(real_p), DIMENSION(nc) :: x,JT
    REAL(real_p), DIMENSION(nc) :: d1,a1,c1,b1,sol1
    REAL(real_p), DIMENSION(mc) :: y,u
    REAL(real_p), DIMENSION(mc) :: d2,a2,c2,b2,sol2
    REAL(real_p), DIMENSION(nc,mc) :: phi, phi_o
    REAL(real_p) :: dx,dy,dxdy,dydx,gdxdy,gdydx,dr2,JB,tran,pi
    REAL(real_p) :: ae,aw,as,an,ao,res,uavg,Pe,res1
    REAL(real_p) :: sum_l,sum_r,sum_b,gam
    REAL(4) :: ta(2),time_start,time_end
    

! Open output file
    OPEN(unit=10, file="Example6_6.out", status = "unknown")
    OPEN(unit=20, file="Example6_6.rsl", status = "unknown")

    pi = 4.0d0*datan(one)
    uavg = 1.0D-3
    gam = kap/(rho*cp)

! x is same as z and y is same as r, height is same as outer radius.
    dx = length/nc; dy = height/mc; dr2 = half*dy
    dydx = dy/dx; dxdy = dx/dy
    gdydx = gam*dydx; gdxdy = gam*dxdy

    x(1) = half*dx
    DO i = 2,nc
      x(i) = x(i-1)+dx
    ENDDO
    y(1) = half*dy
    u(1) = 2.0d0*uavg*(one-(y(1)/height)**2)
    DO j = 2,mc
      y(j) = y(j-1)+dy
      u(j) = 2.0d0*uavg*(one-(y(j)/height)**2)
    ENDDO

! Set initial condition
    phi_o(:,:) = phi_in

! Start of time marching
    DO it = 1, max_times

    phi(:,:) = phi_o(:,:)  ! Set initial guess to previous time step solution

    DO i = 1,nc
      IF(x(i) < 2.0d0*length/3.0d0 .AND. x(i) > length/3.0d0) THEN
        JT(i) = JT0*sin(pi*it*dt/100.0d0)
      ELSE
        JT(i) = zero
      ENDIF
    ENDDO
    
! Start of iteration loop

    time_start = etime(ta)
    
    DO iter = 1,max_iter

! Row-wise sweep
! Bottom Row
      j = 1
      d1(:) = zero; c1(:) = zero; a1(:) = zero; b1(:) = zero
      d1(1) = u(j)*y(j)*dy + 4.0d0*gdydx*y(j) + gdxdy*(y(j)+dr2)
      d1(1) = d1(1) + y(j)*dy*dx/dt
      c1(1) = -4.0d0*gdydx*y(j)/3.0d0
      b1(1) = (u(j)*y(j)*dy+8.0d0*gdydx*y(j)/3.0d0)*phi_in+gdxdy*(y(j)+dr2)*phi(1,j+1)
      b1(1) = b1(1) + y(j)*dy*dx*phi_o(1,j)/dt
      DO i = 2,nc-1
        d1(i) = u(j)*y(j)*dy + 2.0d0*y(j)*gdydx + gdxdy*(y(j)+dr2)
        d1(i) = d1(i) + y(j)*dy*dx/dt
        c1(i) = -gdydx*y(j)
        a1(i-1) = -(u(j)*y(j)*dy+gdydx*y(j))
        b1(i) = gdxdy*(y(j)+dr2)*phi(i,j+1)
        b1(i) = b1(i) + y(j)*dy*dx*phi_o(i,j)/dt
      ENDDO
      d1(nc) = u(j)*y(j)*dy + gdydx*y(j) + gdxdy*(y(j)+dr2)
      d1(nc) = d1(nc) + y(j)*dy*dx/dt
      a1(nc-1) = -(u(j)*y(j)*dy+gdydx*y(j))
      b1(nc) = gdxdy*(y(j)+dr2)*phi(nc,j+1)
      b1(nc) = b1(nc) + y(j)*dy*dx*phi_o(nc,j)/dt

      CALL TRI(nc,a1,d1,c1,b1,sol1)

      phi(:,j) = sol1(:)
     
! Interior Rows

      DO j = 2, mc-1

        d1(:) = zero; c1(:) = zero; a1(:) = zero; b1(:) = zero
        d1(1) = u(j)*y(j)*dy + 4.0d0*y(j)*gdydx + gdxdy*(y(j)+dr2) + gdxdy*(y(j)-dr2)
        d1(1) = d1(1) + y(j)*dy*dx/dt
        c1(1) = -4.0d0*gdydx*y(j)/3.0d0
        b1(1) = (u(j)*y(j)*dy+8.0d0*gdydx*y(j)/3.0d0)*phi_in+ &
                 gdxdy*(y(j)-dr2)*phi(1,j-1)+gdxdy*(y(j)+dr2)*phi(1,j+1)
        b1(1) = b1(1) + y(j)*dy*dx*phi_o(1,j)/dt
        DO i = 2,nc-1
          d1(i) = u(j)*y(j)*dy + 2.0d0*gdydx*y(j) + gdxdy*(y(j)+dr2) + gdxdy*(y(j)-dr2)
          d1(i) = d1(i) + y(j)*dy*dx/dt
          c1(i) = -gdydx*y(j)
          a1(i-1) = -(u(j)*y(j)*dy+gdydx*y(j))
          b1(i) = gdxdy*(y(j)-dr2)*phi(i,j-1)+gdxdy*(y(j)+dr2)*phi(i,j+1)
          b1(i) = b1(i) + y(j)*dy*dx*phi_o(i,j)/dt
        ENDDO
        d1(nc) = u(j)*y(j)*dy + gdydx*y(j) + gdxdy*(y(j)+dr2) + gdxdy*(y(j)-dr2)
        d1(nc) = d1(nc) + y(j)*dy*dx/dt
        a1(nc-1) = -(u(j)*y(j)*dy+gdydx*y(j))
        b1(nc) = gdxdy*(y(j)-dr2)*phi(nc,j-1)+gdxdy*(y(j)+dr2)*phi(nc,j+1)
        b1(nc) = b1(nc) + y(j)*dy*dx*phi_o(nc,j)/dt

        CALL TRI(nc,a1,d1,c1,b1,sol1)

        phi(:,j) = sol1(:)

      ENDDO

! Top Row
      j = mc
      d1(:) = zero; c1(:) = zero; a1(:) = zero; b1(:) = zero
      d1(1) = u(j)*y(j)*dy + 4.0d0*gdydx*y(j) + gdxdy*(y(j)-dr2)
      d1(1) = d1(1) + y(j)*dy*dx/dt
      c1(1) = -4.0d0*gdydx*y(j)/3.0d0
      b1(1) = (u(j)*y(j)*dy+8.0d0*gdydx*y(j)/3.0d0)*phi_in+gdxdy*(y(j)-dr2)*phi(1,j-1)+gam*(y(j)+dr2)*dx*JT(1)/kap
      b1(1) = b1(1) + y(j)*dy*dx*phi_o(1,j)/dt
      DO i = 2,nc-1
        d1(i) = u(j)*y(j)*dy + 2.0d0*y(j)*gdydx + gdxdy*(y(j)-dr2)
        d1(i) = d1(i) + y(j)*dy*dx/dt
        c1(i) = -gdydx*y(j)
        a1(i-1) = -(u(j)*y(j)*dy+gdydx*y(j))
        b1(i) = gdxdy*(y(j)-dr2)*phi(i,j-1)+gam*(y(j)+dr2)*dx*JT(i)/kap
        b1(i) = b1(i) + y(j)*dy*dx*phi_o(i,j)/dt
      ENDDO
      d1(nc) = u(j)*y(j)*dy + gdydx*y(j) + gdxdy*(y(j)-dr2)
      d1(nc) = d1(nc) + y(j)*dy*dx/dt
      a1(nc-1) = -(u(j)*y(j)*dy+gdydx*y(j))
      b1(nc) = gdxdy*(y(j)-dr2)*phi(nc,j-1)+gam*(y(j)+dr2)*dx*JT(nc)/kap
      b1(nc) = b1(nc) + y(j)*dy*dx*phi_o(nc,j)/dt

      CALL TRI(nc,a1,d1,c1,b1,sol1)

      phi(:,j) = sol1(:)

! Residual calculation
      res = zero
! Interior      
      DO i = 2,nc-1
        DO j = 2,mc-1
          ao = (u(j)*y(j)*dy + 2.0d0*gdydx*y(j) + gdxdy*(y(j)+dr2) + gdxdy*(y(j)-dr2))*phi(i,j)
          aw = -(u(j)*y(j)*dy+gdydx*y(j))*phi(i-1,j)
          ae = -gdydx*y(j)*phi(i+1,j)
          as = -gdxdy*(y(j)-dr2)*phi(i,j-1)
          an = -gdxdy*(y(j)+dr2)*phi(i,j+1)
          tran = (phi(i,j)-phi_o(i,j))*y(j)*dy*dx/dt
          res = res + (ao+aw+ae+as+an+tran)**2
        ENDDO
      ENDDO
! Inlet
      i = 1
      DO j = 2,mc-1
        ao = (u(j)*y(j)*dy + 4.0d0*gdydx*y(j) + gdxdy*(y(j)+dr2) + gdxdy*(y(j)-dr2))*phi(i,j)
        aw = -(u(j)*y(j)*dy+8.0d0*gdydx*y(j)/3.0d0)*phi_in
        ae = -4.0d0*gdydx*y(j)*phi(i+1,j)/3.0d0
        as = -gdxdy*(y(j)-dr2)*phi(i,j-1)
        an = -gdxdy*(y(j)+dr2)*phi(i,j+1)
        tran = (phi(i,j)-phi_o(i,j))*y(j)*dy*dx/dt
        res = res + (ao+aw+ae+as+an+tran)**2
      ENDDO
! Outlet
      i = nc
      DO j = 2,mc-1
        ao = (u(j)*y(j)*dy + gdydx*y(j) + gdxdy*(y(j)+dr2) + gdxdy*(y(j)-dr2))*phi(i,j)
        aw = -(u(j)*y(j)*dy+gdydx*y(j))*phi(i-1,j)
        ae = zero
        as = -gdxdy*(y(j)-dr2)*phi(i,j-1)
        an = -gdxdy*(y(j)+dr2)*phi(i,j+1)
        tran = (phi(i,j)-phi_o(i,j))*y(j)*dy*dx/dt
        res = res + (ao+aw+ae+as+an+tran)**2
      ENDDO
! Bottom
      j = 1
      DO i = 2,nc-1
        ao = (u(j)*y(j)*dy + 2.0d0*gdydx*y(j) + gdxdy*(y(j)+dr2))*phi(i,j)
        aw = -(u(j)*y(j)*dy+gdydx*y(j))*phi(i-1,j)
        ae = -gdydx*y(j)*phi(i+1,j)
        as = zero
        an = -gdxdy*(y(j)+dr2)*phi(i,j+1)
        tran = (phi(i,j)-phi_o(i,j))*y(j)*dy*dx/dt
        res = res + (ao+aw+ae+as+an+tran)**2
      ENDDO
! Top
      j = mc
      DO i = 2,nc-1
        ao = (u(j)*y(j)*dy + 2.0d0*gdydx*y(j) + gdxdy*(y(j)-dr2))*phi(i,j)
        aw = -(u(j)*y(j)*dy+gdydx*y(j))*phi(i-1,j)
        ae = -gdydx*y(j)*phi(i+1,j)
        as = -gdxdy*(y(j)-dr2)*phi(i,j-1)
        an = -gam*(y(j)+dr2)*dx*JT(i)/kap
        tran = (phi(i,j)-phi_o(i,j))*y(j)*dy*dx/dt
        res = res + (ao+aw+ae+as+an+tran)**2
      ENDDO

      res = SQRT(MAX(zero,res))
      IF(iter == 1) res1 = res
      WRITE(*,*) iter,res/res1
      WRITE(20,*) iter*2,res/res1

      IF(res/res1 < tol) EXIT

    ENDDO  !iterations

    phi_o(:,:) = phi(:,:)  ! Reset initial condition

    ENDDO  ! Time steps
  
    time_end = etime(ta)
    
! Print results (for last time step)

    WRITE(10,*) 'VARIABLES = "X", "Y", "PHI"'
    WRITE(10,*) 'ZONE I=', nc+2, ', J=', mc+2, ', DATAPACKING=POINT'

!Bottom Row
    WRITE(10,10) zero,zero,phi_in
    DO i = 1,nc
      WRITE(10,10) x(i),zero,(9.0d0*phi(i,1)-phi(i,2))/8.0d0
    ENDDO
    WRITE(10,10) length,zero,(9.0d0*phi(nc,1)-phi(nc,2))/8.0d0

!Interior Rows
    DO j = 1,mc
      WRITE(10,10) zero, y(j), phi_in
      DO i = 1,nc
         WRITE(10,10) x(i),y(j),phi(i,j)
      ENDDO
      WRITE(10,10) length, y(j), phi(nc,j)
    ENDDO

!Top Row
    WRITE(10,10) zero,height,phi_in
    DO i = 1,nc
      WRITE(10,10) x(i),height,(9.0d0*phi(i,mc-1)-phi(i,mc-2)+3.0d0*JT(i)*dy/kap)/8.0d0
    ENDDO
    WRITE(10,10) length,height,(9.0d0*phi(nc,mc-1)-phi(nc,mc-2)+3.0d0*JT(nc)*dy/kap)/8.0d0


    print *, "time taken = ", time_end-time_start

 10 format(3(1x,F13.6))

    CLOSE(unit=10)
    CLOSE(unit=20)

    END Program main

!--------------------------------------------------------------------------

! Tridiagonal solver (uses Thomas Algorithm)

    SUBROUTINE TRI(n,a,d,c,b,x)

    INTEGER(4), INTENT(in) :: n
    INTEGER(4) :: i
    REAL(8) :: xmult
    REAL(8), DIMENSION(n) :: a,b,c,d,x

    DO i = 2,n
      xmult = a(i-1)/d(i-1)
      d(i) = d(i) - xmult*c(i-1)
      b(i) = b(i) - xmult*b(i-1)
    ENDDO

    x(n) = b(n)/d(n)
    DO i = n-1,1,-1
      x(i) = (b(i)-c(i)*x(i+1))/d(i)
    ENDDO

    END SUBROUTINE TRI
