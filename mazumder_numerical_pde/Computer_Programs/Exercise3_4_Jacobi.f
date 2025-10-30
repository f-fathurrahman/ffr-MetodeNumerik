! Start of Main Program
! This is a program to solve a second order linear elliptic PDE
! using the finite-difference method and Jacobi iterations
! Exercise 3.4(a)

    program main

    USE DFPORT

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n = 41, m = 41
    INTEGER(4), PARAMETER :: max_iter = 1000000
    INTEGER(4) :: i,j,iter
    REAL(8), PARAMETER :: one = 1.0D0,two = 2.0d0, &
                             & zero = 0.0d0,half = 0.5d0
    REAL(8), PARAMETER :: length = 1.0D0
    REAL(8) :: tol = 1.0d-6
    REAL(8), DIMENSION(n) :: x,phi_bot,phi_top
    REAL(8), DIMENSION(m) :: y,phi_left,phi_right
    REAL(8), DIMENSION(n,m) :: S,phi,phia,phio
    REAL(8) :: dx,dx2,dy,dy2
    REAL(8) :: ae,aw,as,an,ao,source,res,fde
    REAL(4) :: ta(2),time_start,time_end
    

! Open output file
    OPEN(unit=10, file="Exercise3_4_Jacobi.out", status = "unknown")
    OPEN(unit=20, file="Exercise3_4_Jacobi.rsl", status = "unknown")

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

      phio(:,:) = phi(:,:)

! Jacobi update
      DO i = 2,n-1
        DO j = 2,m-1
          source = -S(i,j)
          phi(i,j) = (source-aw*phio(i-1,j)-ae*phio(i+1,j)-as*phio(i,j-1)-an*phio(i,j+1))/ao
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
      WRITE(*,*) iter,res
      WRITE(20,*) iter,res

      IF(res < tol) EXIT

    ENDDO  !iterations
  
    time_end = etime(ta)
    
! Print results
    WRITE(10,*) 'VARIABLES = "X", "Y", "PHIAN", "PHI", "Error"'
    WRITE(10,*) 'ZONE I=', n, ', J=', m, ', DATAPACKING=POINT'
    DO j = 1,m
      DO i = 1,n
         phia(i,j) = 500.0d0*exp(-50.0*((one-x(i))**2+y(j)**2)) + 100.0d0*x(i)*(one-y(j))
         WRITE(10,10) x(i),y(j),phia(i,j),phi(i,j),phia(i,j)-phi(i,j)
      ENDDO
    ENDDO

    print *, "time taken = ", time_end-time_start

 10 format(5(1x,F13.6))

    CLOSE(unit=10)
    CLOSE(unit=20)

    END Program main
