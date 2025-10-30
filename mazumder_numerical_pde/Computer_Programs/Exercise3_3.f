! Start of Main Program
! Gauss-Seidel Method with Linear Relaxation
! Exercise 3.3

    program main

    USE DFPORT

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n = 81, m = 81
    INTEGER(4), PARAMETER :: max_iter = 1000000
    INTEGER(4) :: i,j,iter
    REAL(8), PARAMETER :: one = 1.0D0,two = 2.0d0, &
                             & zero = 0.0d0,half = 0.5d0
    REAL(8), PARAMETER :: aa = 0.5d0, bb = 0.5d0, cc = 10.0d0, &
                               dd = 1.0d0, ee = 2.0d0
    REAL(8), PARAMETER :: length = 1.0D0
    REAL(8) :: tol = 1.0d-6
    REAL(8), DIMENSION(n) :: x,phi_bot,phi_top
    REAL(8), DIMENSION(m) :: y,phi_left,phi_right
    REAL(8), DIMENSION(n,m) :: S,phi,phian
    REAL(8) :: dx,dx2,dy,dy2,phiold,phinew
    REAL(8) :: ae,aw,as,an,ao,source,res,fde
    REAL(4) :: ta(2),time_start,time_end
    REAL(8), PARAMETER :: omega = 0.5d0
    

! Open output file
    OPEN(unit=10, file="Exercise3_3.out", status = "unknown")
    OPEN(unit=20, file="Exercise3_3.rsl", status = "unknown")

    dx = length/(n-1); dy = length/(m-1)
    dx2 = dx*dx; dy2 = dy*dy
    
    DO i = 1,n
      x(i) = (i-1)*dx
      phi_bot(i) = bb*bb*sinh(-cc*bb)+(x(i)-aa)*(x(i)-aa)*sinh(cc*(x(i)-aa))+dd
      phi_top(i) = (one-bb)*(one-bb)*sinh(cc*(one-bb))+(x(i)-aa)*(x(i)-aa)*sinh(cc*(x(i)-aa))+dd*exp(ee*x(i))
    ENDDO
    DO j = 1,m
      y(j) = (j-1)*dy
      phi_left(j) = aa*aa*sinh(-cc*aa)+(y(j)-bb)*(y(j)-bb)*sinh(cc*(y(j)-bb))+dd
      phi_right(j) = (one-aa)*(one-aa)*sinh(cc*(one-aa))+(y(j)-bb)*(y(j)-bb)*sinh(cc*(y(j)-bb))+dd*exp(ee*y(j))
    ENDDO
    DO i = 1,n
      DO j = 1,m
        S(i,j) = 2.0d0*sinh(cc*(x(i)-aa))+4.0d0*(x(i)-aa)*cc*cosh(cc*(x(i)-aa))+ &
                 (x(i)-aa)*(x(i)-aa)*(cc*cc)*sinh(cc*(x(i)-aa)) + &
                 2.0d0*sinh(cc*(y(j)-bb))+4.0d0*(y(j)-bb)*cc*cosh(cc*(y(j)-bb))+ &
                 (y(j)-bb)*(y(j)-bb)*(cc*cc)*sinh(cc*(y(j)-bb)) + &
                 dd*ee*ee*(y(j)*y(j)+x(i)*x(i))*exp(ee*x(i)*y(j))
      ENDDO
    ENDDO

    phi(:,:) = zero

! Boundary conditions
    DO j = 2,m
      phi(1,j) = phi_left(j)
    ENDDO
    DO j = 2,m
      phi(n,j)= phi_right(j)
    ENDDO
    DO i = 1,n
      phi(i,1) = phi_bot(i)
    ENDDO
    DO i = 2,n-1
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
          phiold = phi(i,j)
          phinew = (source-aw*phi(i-1,j)-ae*phi(i+1,j)-as*phi(i,j-1)-an*phi(i,j+1))/ao
          phi(i,j) = omega*phinew + (one-omega)*phiold
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
         phian(i,j) = ((x(i)-aa)**2)*sinh(cc*(x(i)-aa)) + ((y(j)-bb)**2)*sinh(cc*(y(j)-bb))+dd*exp(ee*x(i)*y(j))
         WRITE(10,10) x(i),y(j),phian(i,j),phi(i,j),phian(i,j)-phi(i,j)
      ENDDO
    ENDDO

    print *, "time taken = ", time_end-time_start

 10 format(5(1x,F13.6))

    CLOSE(unit=10)
    CLOSE(unit=20)

    END Program main
