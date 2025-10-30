! Start of Main Program
! This is a program to solve a second order linear elliptic PDE
! using the finite-difference method and Gauss-Seidel iterations
! Exercise 6.5(a) [same as Exercise 3.4(b)]

    program main

    USE DFPORT

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n = 81, m = 81
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
    REAL(8) :: ae,aw,as,an,ao,source,res,fde,res1
    REAL(8) :: sumf_b,sumf_t,sumf_l,sumf_r,imb,prod
    REAL(4) :: ta(2),time_start,time_end
    

! Open output file
    OPEN(unit=10, file="Exercise6_5_FDM.out", status = "unknown")
    OPEN(unit=20, file="Exercise6_5_FDM.rsl", status = "unknown")
    OPEN(unit=30, file="Exercise6_5_FDM.flx", status = "unknown")

    dx = length/(n-1); dy = length/(m-1)
    dx2 = dx*dx; dy2 = dy*dy
    
    DO i = 1,n
      x(i) = (i-1)*dx
      phi_bot(i) = 100.0d0*x(i)+500.0d0*exp(-50.0*(one-x(i))**2) !bottom
      phi_top(i) = 500.0d0*exp(-50.0*((one-x(i))**2+one)) !top
    ENDDO
    DO j = 1,m
      y(j) = (j-1)*dy
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

! Boundary conditions
    DO j = 2,m-1
      phi(1,j) = phi_left(j)
    ENDDO
    DO j = 2,m-1
      phi(n,j)= phi_right(j)
    ENDDO
    DO i = 1,n
      phi(i,1) = phi_bot(i)
    ENDDO
    DO i = 1,n
      phi(i,m) = phi_top(i)
    ENDDO


! Link Coefficients
    aw = -one/dx2
    ae = -one/dx2
    as = -one/dy2
    an = -one/dy2
    ao = two/dx2+two/dy2

! Start of iteration loop

    time_start = etime(ta)
    

    DO iter = 1,max_iter

! Gauss-Seidel update
      DO i = 2,n-1
        DO j = 2,m-1
          source = -S(i,j)
          phi(i,j) = (source-aw*phi(i-1,j)-ae*phi(i+1,j)-as*phi(i,j-1)-an*phi(i,j+1))/ao
        ENDDO
      ENDDO

! Residual calculation
      res = zero      
      DO i = 2,n-1
        DO j = 2,m-1
          source = -S(i,j)
          fde = source-aw*phi(i-1,j)-ae*phi(i+1,j)-as*phi(i,j-1)-an*phi(i,j+1)-ao*phi(i,j)
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
    WRITE(10,*) 'ZONE I=', n, ', J=', m, ', DATAPACKING=POINT'
    DO j = 1,m
      DO i = 1,n
         phia(i,j) = 500.0d0*exp(-50.0*((one-x(i))**2+y(j)**2)) + 100.0d0*x(i)*(one-y(j))
         WRITE(10,10) x(i),y(j),phia(i,j),phi(i,j),phia(i,j)-phi(i,j)
      ENDDO
    ENDDO
! Compute and Print results for Fluxes
    sumf_b = zero
    sumf_t = zero
    WRITE(30,*) " x(i )", " Flux Bottom ", " Flux Top "
    DO i = 1,n
      flux_b(i) = (4.0d0*phi(i,2)-phi(i,3)-3.0d0*phi(i,1))/(2.0d0*dy)
      flux_t(i) = -(4.0d0*phi(i,m-1)-phi(i,m-2)-3.0d0*phi(i,m))/(2.0d0*dy)
      IF(i == 1 .OR. i == n)THEN
        sumf_b = sumf_b + flux_b(i)*dx/2.0d0
        sumf_t = sumf_t + flux_t(i)*dx/2.0d0
      ELSE
        sumf_b = sumf_b + flux_b(i)*dx
        sumf_t = sumf_t + flux_t(i)*dx
      ENDIF
      WRITE(30,*) x(i),flux_b(i),flux_t(i)
      WRITE(*,*) x(i),flux_b(i),flux_t(i)
    ENDDO
    sumf_l = zero
    sumf_r = zero
    WRITE(30,*) " y(i) ", " Flux Left ", " Flux Right "
    DO j = 1,m
      flux_l(j) = (4.0d0*phi(2,j)-phi(3,j)-3.0d0*phi(1,j))/(2.0d0*dx)
      flux_r(j) = -(4.0d0*phi(n-1,j)-phi(n-2,j)-3.0d0*phi(n,j))/(2.0d0*dx)
      IF(j == 1 .OR. j == m)THEN
        sumf_l = sumf_l + flux_l(j)*dy/2.0d0
        sumf_r = sumf_r + flux_r(j)*dy/2.0d0
      ELSE
        sumf_l = sumf_l + flux_l(j)*dy
        sumf_r = sumf_r + flux_r(j)*dy
      ENDIF
      WRITE(30,*) y(j),flux_l(j),flux_r(j)
      WRITE(*,*) y(j),flux_l(j),flux_r(j)
    ENDDO

    prod = zero
    prod = prod + S(1,1)*dx*dy/4.0d0
    DO i = 2,n-1
      prod = prod + S(i,1)*dx*dy/2.0d0
    ENDDO
    prod = prod + S(n,1)*dx*dy/4.0d0
    DO j = 2,m-1
      prod = prod + S(n,j)*dx*dy/2.0d0
    ENDDO
    prod = prod + S(n,m)*dx*dy/4.0d0
    DO i = 2,n-1
      prod = prod + S(i,m)*dx*dy/2.0d0
    ENDDO
    prod = prod + S(1,m)*dx*dy/4.0d0
    DO j = 2,m-1
      prod = prod + S(1,j)*dx*dy/2.0d0
    ENDDO
    DO i = 2,n-1
      DO j = 2,m-1
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
