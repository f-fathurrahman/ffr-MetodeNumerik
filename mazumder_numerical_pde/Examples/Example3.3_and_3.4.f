! Start of Main Program
! Example 3.3 and 3.4
! ADI Method

    program main

    USE DFPORT

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n = 41, m = 81
    INTEGER(4), PARAMETER :: max_iter = 1000000
    INTEGER(4) :: i,j,iter
    REAL(8), PARAMETER :: one = 1.0D0,two = 2.0d0, &
                             & zero = 0.0d0,half = 0.5d0
    REAL(8), PARAMETER :: aa = 0.5d0, bb = 0.5d0, cc = 10.0d0, &
                               dd = 1.0d0, ee = 2.0d0
    REAL(8), PARAMETER :: length = 1.0D0
    REAL(8) :: tol = 1.0d-6
    REAL(8), DIMENSION(n) :: x,phi_bot,phi_top
    REAL(8), DIMENSION(n) :: d1,a1,c1,b1,sol1
    REAL(8), DIMENSION(m) :: y,phi_left,phi_right
    REAL(8), DIMENSION(m) :: d2,a2,c2,b2,sol2
    REAL(8), DIMENSION(n,m) :: S,phi,phian
    REAL(8) :: dx,dx2,dy,dy2
    REAL(8) :: ae,aw,as,an,ao,source,res,fde
    REAL(4) :: ta(2),time_start,time_end
    

! Open output file
    OPEN(unit=10, file="Ch3_ADI.out", status = "unknown")
    OPEN(unit=20, file="Ch3_ADI.rsl", status = "unknown")

    dx = length/(n-1); dy = length/(m-1)
    dx2 = dx*dx; dy2 = dy*dy
    
    DO i = 1,n
      x(i) = (i-1)*dx
      phi_bot(i) = bb*bb*sinh(-cc*bb)+(x(i)-aa)*(x(i)-aa)*sinh(cc*(x(i)-aa))+dd
      phi_top(i) = (one-bb)*(one-bb)*sinh(cc*(one-bb))+(x(i)-aa)* &
                   (x(i)-aa)*sinh(cc*(x(i)-aa))+dd*exp(ee*x(i))
    ENDDO
    DO j = 1,m
      y(j) = (j-1)*dy
      phi_left(j) = aa*aa*sinh(-cc*aa)+(y(j)-bb)*(y(j)-bb)*sinh(cc*(y(j)-bb))+dd
      phi_right(j) = (one-aa)*(one-aa)*sinh(cc*(one-aa))+ &
                     (y(j)-bb)*(y(j)-bb)*sinh(cc*(y(j)-bb))+dd*exp(ee*y(j))
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

! Row-wise sweep
      DO j = 2,m-1
        
        d1(:) = zero; c1(:) = zero; a1(:) = zero; b1(:) = zero
        d1(1) = one
        b1(1) = phi_left(j)
        DO i = 2,n-1
          d1(i) = ao
          c1(i) = ae
          a1(i-1) = aw
          b1(i) = -S(i,j)-an*phi(i,j+1)-as*phi(i,j-1)
        ENDDO
        d1(n) = one
        b1(n) = phi_right(j)

        CALL TRI(n,a1,d1,c1,b1,sol1)

        phi(:,j) = sol1(:)

      ENDDO

! Column-wise sweep
      DO i = 2,n-1
        
        d2(:) = zero; c2(:) = zero; a2(:) = zero; b2(:) = zero
        d2(1) = one
        b2(1) = phi_bot(i)
        DO j = 2,m-1
          d2(j) = ao
          c2(j) = an
          a2(j-1) = as
          b2(j) = -S(i,j)-ae*phi(i+1,j)-aw*phi(i-1,j)
        ENDDO
        d2(m) = one
        b2(m) = phi_top(i)

        CALL TRI(m,a2,d2,c2,b2,sol2)

        phi(i,:) = sol2(:)

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
      WRITE(*,*) iter*2,res
      WRITE(20,*) iter*2,res

      IF(res < tol) EXIT

    ENDDO  !iterations
  
    time_end = etime(ta)
    
! Print results
    WRITE(10,*) 'VARIABLES = "X", "Y", "PHIAN", "PHI", "Error"'
    WRITE(10,*) 'ZONE I=', n, ', J=', m, ', DATAPACKING=POINT'
    DO j = 1,m
      DO i = 1,n
         phian(i,j) = ((x(i)-aa)**2)*sinh(cc*(x(i)-aa)) + &
                      ((y(j)-bb)**2)*sinh(cc*(y(j)-bb))+dd*exp(ee*x(i)*y(j))
         WRITE(10,10) x(i),y(j),phian(i,j),phi(i,j),phian(i,j)-phi(i,j)
      ENDDO
    ENDDO

    print *, "time taken = ", time_end-time_start

 10 format(5(1x,F13.6))

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
