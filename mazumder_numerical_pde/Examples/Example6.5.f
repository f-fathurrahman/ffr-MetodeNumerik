! Start of Main Program
! Example 6.5
! ADI Method

    program main

    USE DFPORT

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: nc = 200, mc = 40
    INTEGER(4), PARAMETER :: max_iter = 1000000
    INTEGER(4) :: i,j,iter
    REAL(8), PARAMETER :: one = 1.0D0,two = 2.0d0, &
                             & zero = 0.0d0,half = 0.5d0
    REAL(8), PARAMETER :: length = 0.5D0, height = 0.1d0, rho = 1000.0d0, &
                               JB = 10.0d0, kap = 0.613d0, gam = 1.45D-7, &
                               phi_in = 20.0d0
    REAL(8) :: tol = 1.0d-6
    REAL(8), DIMENSION(nc) :: x,phi_bot,phi_top
    REAL(8), DIMENSION(nc) :: d1,a1,c1,b1,sol1
    REAL(8), DIMENSION(mc) :: y,u
    REAL(8), DIMENSION(mc) :: d2,a2,c2,b2,sol2
    REAL(8), DIMENSION(nc,mc) :: phi
    REAL(8) :: dx,dy,dxdy,dydx,gdxdy,gdydx
    REAL(8) :: ae,aw,as,an,ao,res,uavg,Pe,res1
    REAL(8) :: sum_l,sum_r,sum_b,cp
    REAL(4) :: ta(2),time_start,time_end
    

! Open output file
    OPEN(unit=10, file="Example6_5.out", status = "unknown")
    OPEN(unit=20, file="Example6_5.rsl", status = "unknown")

    dx = length/nc; dy = height/mc
    dydx = dy/dx; dxdy = dx/dy
    gdydx = gam*dydx; gdxdy = gam*dxdy

    Pe = 200.0d0
    uavg = Pe*gam/(rho*height)
    cp = kap/(rho*gam)
    
    x(1) = half*dx
    DO i = 2,nc
      x(i) = x(i-1)+dx
    ENDDO
    y(1) = half*dy
    u(1) = 6.0d0*uavg*((y(1)/height)-(y(1)/height)**2)
    DO j = 2,mc
      y(j) = y(j-1)+dy
      u(j) = 6.0d0*uavg*((y(j)/height)-(y(j)/height)**2)
    ENDDO

    phi(:,:) = phi_in

! Start of iteration loop

    time_start = etime(ta)
    
    DO iter = 1,max_iter

! Row-wise sweep
! Bottom Row
      j = 1
      d1(:) = zero; c1(:) = zero; a1(:) = zero; b1(:) = zero
      d1(1) = u(j)*dy + 4.0d0*gdydx + gdxdy
      c1(1) = -4.0d0*gdydx/3.0d0
      b1(1) = (u(j)*dy+8.0d0*gdydx/3.0d0)*phi_in+gam*JB*dx/kap+gdxdy*phi(1,j+1)
      DO i = 2,nc-1
        d1(i) = u(j)*dy + 2.0d0*gdydx + gdxdy
        c1(i) = -gdydx
        a1(i-1) = -(u(j)*dy+gdydx)
        b1(i) = gam*JB*dx/kap+gdxdy*phi(i,j+1)
      ENDDO
      d1(nc) = u(j)*dy + gdydx + gdxdy
      a1(nc-1) = -(u(j)*dy+gdydx)
      b1(nc) = gam*JB*dx/kap+gdxdy*phi(nc,j+1)

      CALL TRI(nc,a1,d1,c1,b1,sol1)

      phi(:,j) = sol1(:)
     
! Interior Rows

      DO j = 2, mc-1

        d1(:) = zero; c1(:) = zero; a1(:) = zero; b1(:) = zero
        d1(1) = u(j)*dy + 4.0d0*gdydx + 2.0d0*gdxdy
        c1(1) = -4.0d0*gdydx/3.0d0
        b1(1) = (u(j)*dy+8.0d0*gdydx/3.0d0)*phi_in+gdxdy*phi(1,j-1)+gdxdy*phi(1,j+1)
        DO i = 2,nc-1
          d1(i) = u(j)*dy + 2.0d0*gdydx + 2.0d0*gdxdy
          c1(i) = -gdydx
          a1(i-1) = -(u(j)*dy+gdydx)
          b1(i) = gdxdy*phi(i,j-1)+gdxdy*phi(i,j+1)
        ENDDO
        d1(nc) = u(j)*dy + gdydx + 2.0d0*gdxdy
        a1(nc-1) = -(u(j)*dy+gdydx)
        b1(nc) = gdxdy*phi(nc,j-1)+gdxdy*phi(nc,j+1)

        CALL TRI(nc,a1,d1,c1,b1,sol1)

        phi(:,j) = sol1(:)

      ENDDO

! Top Row
      j = mc
      d1(:) = zero; c1(:) = zero; a1(:) = zero; b1(:) = zero
      d1(1) = u(j)*dy + 4.0d0*gdydx + gdxdy
      c1(1) = -4.0d0*gdydx/3.0d0
      b1(1) = (u(j)*dy+8.0d0*gdydx/3.0d0)*phi_in+gdxdy*phi(1,j-1)
      DO i = 2,nc-1
        d1(i) = u(j)*dy + 2.0d0*gdydx + gdxdy
        c1(i) = -gdydx
        a1(i-1) = -(u(j)*dy+gdydx)
        b1(i) = gdxdy*phi(i,j-1)
      ENDDO
      d1(nc) = u(j)*dy + gdydx + gdxdy
      a1(nc-1) = -(u(j)*dy+gdydx)
      b1(nc) = gdxdy*phi(nc,j-1)

      CALL TRI(nc,a1,d1,c1,b1,sol1)

      phi(:,j) = sol1(:)

! Column-wise sweep
! Left Column
      i = 1
      d2(:) = zero; c2(:) = zero; a2(:) = zero; b2(:) = zero
      d2(1) = u(1)*dy + 4.0d0*gdydx + gdxdy
      c2(1) = -gdxdy
      b2(1) = (u(1)*dy+8.0d0*gdydx/3.0d0)*phi_in+gam*JB*dx/kap+4.0d0*gdydx*phi(i+1,1)/3.0d0
      DO j = 2,mc-1
        d2(j) = u(j)*dy + 4.0d0*gdydx + 2.0d0*gdxdy
        c2(j) = -gdxdy
        a2(j-1) = -gdxdy
        b2(j) = 4.0d0*gdydx*phi(i+1,j)/3.0d0+(u(j)*dy+8.0d0*gdydx/3.0d0)*phi_in
      ENDDO
      d2(mc) = u(mc)*dy + 4.0d0*gdydx + gdxdy
      a2(mc-1) = -gdxdy
      b2(mc) = 4.0d0*gdydx*phi(i+1,mc)/3.0d0+(u(mc)*dy+8.0d0*gdydx/3.0d0)*phi_in

      CALL TRI(mc,a2,d2,c2,b2,sol2)

      phi(i,:) = sol2(:)
     
! Interior Columns

      DO i = 2, nc-1

        d2(:) = zero; c2(:) = zero; a2(:) = zero; b2(:) = zero
        d2(1) = u(1)*dy + 2.0d0*gdydx + gdxdy
        c2(1) = -gdxdy
        b2(1) = (u(1)*dy+gdydx)*phi(i-1,1)+gdydx*phi(i+1,1)+gam*JB*dx/kap
        DO j = 2,mc-1
          d2(j) = u(j)*dy + 2.0d0*gdydx + 2.0d0*gdxdy
          c2(j) = -gdxdy
          a2(j-1) = -gdxdy
          b2(j) = gdydx*phi(i+1,j)+(u(j)*dy+gdydx)*phi(i-1,j)
        ENDDO
        d2(mc) = u(mc)*dy + gdxdy + 2.0d0*gdydx
        a2(mc-1) = -gdxdy
        b2(mc) = gdydx*phi(i+1,mc)+(u(mc)*dy+gdydx)*phi(i-1,mc)

        CALL TRI(mc,a2,d2,c2,b2,sol2)

        phi(i,:) = sol2(:)

      ENDDO

! Right Column
      i = nc
      d2(:) = zero; c2(:) = zero; a2(:) = zero; b2(:) = zero
      d2(1) = u(1)*dy + gdydx + gdxdy
      c2(1) = -gdxdy
      b2(1) = (u(1)*dy+gdydx)*phi(i-1,1)+gam*JB*dx/kap
      DO j = 2,mc-1
        d2(j) = u(j)*dy + gdydx + 2.0d0*gdxdy
        c2(j) = -gdxdy
        a2(j-1) = -gdxdy
        b2(j) = (u(j)*dy+gdydx)*phi(i-1,j)
      ENDDO
      d2(mc) = u(mc)*dy + gdydx + gdxdy
      a2(mc-1) = -gdxdy
      b2(mc) = (u(mc)*dy+gdydx)*phi(i-1,mc)

      CALL TRI(mc,a2,d2,c2,b2,sol2)

      phi(i,:) = sol2(:)


! Residual calculation
      res = zero
! Interior      
      DO i = 2,nc-1
        DO j = 2,mc-1
          ao = (u(j)*dy + 2.0d0*gdydx + 2.0d0*gdxdy)*phi(i,j)
          aw = -(u(j)*dy+gdydx)*phi(i-1,j)
          ae = -gdydx*phi(i+1,j)
          as = -gdxdy*phi(i,j-1)
          an = -gdxdy*phi(i,j+1)
          res = res + (ao+aw+ae+as+an)**2
        ENDDO
      ENDDO
! Inlet
      i = 1
      DO j = 2,mc-1
        ao = (u(j)*dy + 4.0d0*gdydx + 2.0d0*gdxdy)*phi(i,j)
        aw = -(u(j)*dy+8.0d0*gdydx/3.0d0)*phi_in
        ae = -4.0d0*gdydx*phi(i+1,j)/3.0d0
        as = -gdxdy*phi(i,j-1)
        an = -gdxdy*phi(i,j+1)
        res = res + (ao+aw+ae+as+an)**2
      ENDDO
! Outlet
      i = nc
      DO j = 2,mc-1
        ao = (u(j)*dy + gdydx + 2.0d0*gdxdy)*phi(i,j)
        aw = -(u(j)*dy+gdydx)*phi(i-1,j)
        ae = zero
        as = -gdxdy*phi(i,j-1)
        an = -gdxdy*phi(i,j+1)
        res = res + (ao+aw+ae+as+an)**2
      ENDDO
! Bottom
      j = 1
      DO i = 2,nc-1
        ao = (u(j)*dy + 2.0d0*gdydx + gdxdy)*phi(i,j)
        aw = -(u(j)*dy+gdydx)*phi(i-1,j)
        ae = -gdydx*phi(i+1,j)
        as = -gam*JB*dx/kap
        an = -gdxdy*phi(i,j+1)
        res = res + (ao+aw+ae+as+an)**2
      ENDDO
! Top
      j = mc
      DO i = 2,nc-1
        ao = (u(j)*dy + 2.0d0*gdydx + gdxdy)*phi(i,j)
        aw = -(u(j)*dy+gdydx)*phi(i-1,j)
        ae = -gdydx*phi(i+1,j)
        as = -gdxdy*phi(i,j-1)
        an = zero
        res = res + (ao+aw+ae+as+an)**2
      ENDDO

      res = SQRT(MAX(zero,res))
      IF(iter == 1) res1 = res
      WRITE(*,*) iter*2,res/res1
      WRITE(20,*) iter*2,res/res1

      IF(res/res1 < tol) EXIT

    ENDDO  !iterations
  
    time_end = etime(ta)
    
! Print results

    WRITE(10,*) 'VARIABLES = "X", "Y", "PHI"'
    WRITE(10,*) 'ZONE I=', nc+2, ', J=', mc+2, ', DATAPACKING=POINT'

!Bottom Row
    WRITE(10,10) zero,zero,phi_in
    DO i = 1,nc
      WRITE(10,10) x(i),zero,(9.0d0*phi(i,1)-phi(i,2)+3.0d0*JB*dy/kap)/8.0d0
    ENDDO
    WRITE(10,10) length,zero,(9.0d0*phi(nc,1)-phi(nc,2)+3.0d0*JB*dy/kap)/8.0d0

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
      WRITE(10,10) x(i),height,(9.0d0*phi(i,mc-1)-phi(i,mc-2))/8.0d0
    ENDDO
    WRITE(10,10) length,height,(9.0d0*phi(nc,mc-1)-phi(nc,mc-2))/8.0d0

! Compute Fluxes
    sum_l = zero  !Inlet
    DO j = 1,mc
      sum_l = sum_l + (rho*u(j)*cp*phi_in - kap*(9.0d0*phi(1,j)-phi(2,j)-8.0d0*phi_in)/(3.0d0*dx))*dy
    ENDDO
    sum_r = zero  !Outlet
    DO j = 1,mc
      sum_r = sum_r + rho*u(j)*cp*phi(nc,j)*dy
    ENDDO
    sum_b = zero  !Bottom
    DO i = 1,nc
      sum_b = sum_b + JB*dx
    ENDDO
    Print *, sum_l,sum_r,sum_b, (sum_r-sum_l-sum_b)

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
