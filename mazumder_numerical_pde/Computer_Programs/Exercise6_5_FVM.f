! Start of Main Program
! This is a program to solve a second order linear elliptic PDE
! using the finite-volume method and Gauss-Seidel iterations
! Exercise 6.5(b)

    program main

    USE DFPORT

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n = 80, m = 80
    INTEGER(4), PARAMETER :: max_iter = 1000000
    INTEGER(4) :: i,j,iter
    REAL(8), PARAMETER :: one = 1.0D0,two = 2.0d0, &
                             & zero = 0.0d0,half = 0.5d0
    REAL(8), PARAMETER :: length = 1.0D0
    REAL(8) :: tol = 1.0d-10
    REAL(8), DIMENSION(n) :: x,phi_bot,phi_top,flux_t,flux_b
    REAL(8), DIMENSION(m) :: y,phi_left,phi_right,flux_l,flux_r
    REAL(8), DIMENSION(n,m) :: S,phi,phia
    REAL(8) :: dx,dx2,dy,dy2
    REAL(8), DIMENSION(n,m) :: ao,ae,aw,as,an,qq
    REAL(8) :: res,fde,res1,dxdy,dydx
    REAL(8) :: sumf_b,sumf_t,sumf_l,sumf_r,imb,prod
    REAL(8) :: phi_l_t,phi_r_t,phi_l_b,phi_r_b
    REAL(4) :: ta(2),time_start,time_end
    

! Open output file
    OPEN(unit=10, file="Exercise6_5_FVM.out", status = "unknown")
    OPEN(unit=20, file="Exercise6_5_FVM.rsl", status = "unknown")
    OPEN(unit=30, file="Exercise6_5_FVM.flx", status = "unknown")

    dx = length/(n); dy = length/(m)
    dxdy = dx/dy; dydx = dy/dx
    x(1) = 0.5d0*dx
    y(1) = 0.5d0*dy
    DO i = 2,n
      x(i) = x(i-1)+dx
    ENDDO
    DO j = 2,m
      y(j) = y(j-1)+dy
    ENDDO
    DO i = 1,n
      phi_bot(i) = 100.0d0*x(i)+500.0d0*exp(-50.0*(one-x(i))**2) !bottom
      phi_top(i) = 500.0d0*exp(-50.0*((one-x(i))**2+one)) !top
    ENDDO
    DO j = 1,m
      phi_left(j) = 500.0d0*exp(-50.0*(one+y(j)**2)) !left
      phi_right(j) = 100*(one-y(j))+500.0d0*exp(-50.0*y(j)**2) !right
    ENDDO
    DO i = 1,n
      DO j = 1,m
        S(i,j) = 5.0D4*(100.0d0*((one-x(i))**2+y(j)**2)-two)* &
                exp(-50.0d0*((one-x(i))**2+y(j)**2))
      ENDDO
    ENDDO

    phi(:,:) = zero

! Link Coefficients and Sources
! Interior cells
    DO i = 2,n-1
      DO j = 2,m-1
        ao(i,j) = 2.0d0*(dxdy+dydx)
        ae(i,j) = -dydx
        aw(i,j) = -dydx
        an(i,j) = -dxdy
        as(i,j) = -dxdy
        qq(i,j) = -S(i,j)*dx*dy
      ENDDO
    ENDDO
! Left boundary
    i = 1
    DO j = 2,m-1
      ao(i,j) = 2.0d0*dxdy+4.0d0*dydx
      ae(i,j) = -4.0d0*dydx/3.0d0
      aw(i,j) = zero
      an(i,j) = -dxdy
      as(i,j) = -dxdy
      qq(i,j) = -S(i,j)*dx*dy + 8.0d0*dydx*phi_left(j)/3.0d0
    ENDDO
! Right boundary
    i = n
    DO j = 2,m-1
      ao(i,j) = 2.0d0*dxdy+4.0d0*dydx
      ae(i,j) = zero
      aw(i,j) = -4.0d0*dydx/3.0d0
      an(i,j) = -dxdy
      as(i,j) = -dxdy
      qq(i,j) = -S(i,j)*dx*dy + 8.0d0*dydx*phi_right(j)/3.0d0
    ENDDO
! Bottom boundary
    j = 1
    DO i = 2,n-1
      ao(i,j) = 4.0d0*dxdy+2.0d0*dydx
      ae(i,j) = -dydx
      aw(i,j) = -dydx
      an(i,j) = -4.0d0*dxdy/3.0d0
      as(i,j) = zero
      qq(i,j) = -S(i,j)*dx*dy + 8.0d0*dxdy*phi_bot(i)/3.0d0
    ENDDO
! Top boundary
    j = m
    DO i = 2,n-1
      ao(i,j) = 4.0d0*dxdy+2.0d0*dydx
      ae(i,j) = -dydx
      aw(i,j) = -dydx
      an(i,j) = zero
      as(i,j) = -4.0d0*dxdy/3.0d0
      qq(i,j) = -S(i,j)*dx*dy + 8.0d0*dxdy*phi_top(i)/3.0d0
    ENDDO
! Left bottom corner
    i = 1; j = 1
    ao(i,j) = 4.0d0*dxdy+4.0d0*dydx
    ae(i,j) = -4.0d0*dydx/3.0d0
    aw(i,j) = zero
    an(i,j) = -4.0d0*dxdy/3.0d0
    as(i,j) = zero
    qq(i,j) = -S(i,j)*dx*dy + 8.0d0*dydx*phi_left(j)/3.0d0 + 8.0d0*dxdy*phi_bot(i)/3.0d0
! Right bottom corner
    i = n; j = 1
    ao(i,j) = 4.0d0*dxdy+4.0d0*dydx
    ae(i,j) = zero
    aw(i,j) = -4.0d0*dydx/3.0d0
    an(i,j) = -4.0d0*dxdy/3.0d0
    as(i,j) = zero
    qq(i,j) = -S(i,j)*dx*dy + 8.0d0*dydx*phi_right(j)/3.0d0 + 8.0d0*dxdy*phi_bot(i)/3.0d0
! Left top corner
    i = 1; j = m
    ao(i,j) = 4.0d0*dxdy+4.0d0*dydx
    ae(i,j) = -4.0d0*dydx/3.0d0
    aw(i,j) = zero
    an(i,j) = zero
    as(i,j) = -4.0d0*dxdy/3.0d0
    qq(i,j) = -S(i,j)*dx*dy + 8.0d0*dydx*phi_left(j)/3.0d0 + 8.0d0*dxdy*phi_top(i)/3.0d0
! Right top corner
    i = n; j = m
    ao(i,j) = 4.0d0*dxdy+4.0d0*dydx
    ae(i,j) = zero
    aw(i,j) = -4.0d0*dydx/3.0d0
    an(i,j) = zero
    as(i,j) = -4.0d0*dxdy/3.0d0
    qq(i,j) = -S(i,j)*dx*dy + 8.0d0*dydx*phi_right(j)/3.0d0 + 8.0d0*dxdy*phi_top(i)/3.0d0

    
! Start of iteration loop

    time_start = etime(ta)
    

    DO iter = 1,max_iter

! Gauss-Seidel update
      DO i = 1,n
        DO j = 1,m
          phi(i,j) = (qq(i,j)-aw(i,j)*phi(i-1,j)-ae(i,j)*phi(i+1,j)- &
                      as(i,j)*phi(i,j-1)-an(i,j)*phi(i,j+1))/ao(i,j)
        ENDDO
      ENDDO

! Residual calculation
      res = zero      
      DO i = 1,n
        DO j = 1,m
          fde = qq(i,j)-aw(i,j)*phi(i-1,j)-ae(i,j)*phi(i+1,j)-as(i,j)*phi(i,j-1)- &
                an(i,j)*phi(i,j+1)-ao(i,j)*phi(i,j)
          res = res + fde*fde
        ENDDO
      ENDDO

      res = SQRT(MAX(zero,res))
      IF(iter == 1) res1 = res
      WRITE(*,*) iter,res/res1
      WRITE(20,*) iter,res/res1

      IF(res/res1 < tol) EXIT

    ENDDO  !iterations
  
    time_end = etime(ta)
    
! Print results for Contour
    WRITE(10,*) 'VARIABLES = "X", "Y", "PHIAN", "PHI", "Error"'
    WRITE(10,*) 'ZONE I=', n+2, ', J=', m+2, ', DATAPACKING=POINT'
    phi_l_b = 500.0d0*exp(-50.0*((one-zero)**2+zero**2)) + 100.0d0*zero*(one-zero)
    WRITE(10,10) zero,zero,phi_l_b,phi_l_b,zero
    DO i = 1,n
      phia(i,j) = 500.0d0*exp(-50.0*((one-x(i))**2+zero**2)) + 100.0d0*x(i)*(one-zero)
      WRITE(10,10) x(i),zero,phia(i,j),phi_bot(i),zero
    ENDDO
    phi_r_b = 500.0d0*exp(-50.0*((one-length)**2+zero**2)) + 100.0d0*length*(one-zero)
    WRITE(10,10) length,zero,phi_r_b,phi_r_b,zero
    DO j = 1,m
      WRITE(10,10) zero,y(j),phi_left(j),phi_left(j),zero
      DO i = 1,n
         phia(i,j) = 500.0d0*exp(-50.0*((one-x(i))**2+y(j)**2)) + 100.0d0*x(i)*(one-y(j))
         WRITE(10,10) x(i),y(j),phia(i,j),phi(i,j),phia(i,j)-phi(i,j)
      ENDDO
      WRITE(10,10) length,y(j),phi_right(j),phi_right(j),zero
    ENDDO
    phi_l_t = 500.0d0*exp(-50.0*((one-zero)**2+length**2)) + 100.0d0*zero*(one-length)
    WRITE(10,10) zero,length,phi_l_t,phi_l_t,zero
    DO i = 1,n
      phia(i,j) = 500.0d0*exp(-50.0*((one-x(i))**2+length**2)) + 100.0d0*x(i)*(one-length)
      WRITE(10,10) x(i),length,phia(i,j),phi_top(i),zero
    ENDDO
    phi_r_t = 500.0d0*exp(-50.0*((one-length)**2+length**2)) + 100.0d0*length*(one-length)
    WRITE(10,10) length,length,phi_r_t,phi_r_t,zero

! Compute and Print results for Fluxes
    sumf_b = zero
    sumf_t = zero
    WRITE(30,*) " x(i )", " Flux Bottom ", " Flux Top "
    DO i = 1,n
      flux_b(i) = (9.0d0*phi(i,1)-phi(i,2)-8.0d0*phi_bot(i))/(3.0d0*dy)
      flux_t(i) = -(9.0d0*phi(i,m)-phi(i,m-1)-8.0d0*phi_top(i))/(3.0d0*dy)
      sumf_b = sumf_b + flux_b(i)*dx
      sumf_t = sumf_t + flux_t(i)*dx
      WRITE(30,*) x(i),flux_b(i),flux_t(i)
      WRITE(*,*) x(i),flux_b(i),flux_t(i)
    ENDDO
    sumf_l = zero
    sumf_r = zero
    WRITE(30,*) " y(i) ", " Flux Left ", " Flux Right "
    DO j = 1,m
      flux_l(j) = (9.0d0*phi(1,j)-phi(2,j)-8.0d0*phi_left(j))/(3.0d0*dx)
      flux_r(j) = -(9.0d0*phi(n,j)-phi(n-1,j)-8.0d0*phi_right(j))/(3.0d0*dx)
      sumf_l = sumf_l + flux_l(j)*dy
      sumf_r = sumf_r + flux_r(j)*dy
      WRITE(30,*) y(j),flux_l(j),flux_r(j)
      WRITE(*,*) y(j),flux_l(j),flux_r(j)
    ENDDO

    prod = zero
    DO i = 1,n
      DO j = 1,m
        prod = prod + S(i,j)*dx*dy
      ENDDO
    ENDDO

    imb = (sumf_r-sumf_l) + (sumf_t-sumf_b)
    WRITE(30,*) " Left ", " Right ", " Bottom ", " Top ", " Production ", " Imbalance "
    WRITE(30,*) sumf_l,sumf_r,sumf_b,sumf_t,prod,imb
    WRITE(*,*) sumf_l,sumf_r,sumf_b,sumf_t,prod,imb



    print *, "time taken = ", time_end-time_start

 10 format(5(1x,F13.6))

    CLOSE(unit=10)
    CLOSE(unit=20)
    CLOSE(unit=30)

    END Program main
