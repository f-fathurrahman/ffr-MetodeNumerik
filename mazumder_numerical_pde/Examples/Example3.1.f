! Start of Main Program
! Example 3.1
program main
  !USE DFPORT

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n = 20, m = 20, neq = n*m
    INTEGER(4) :: i,j,k
    REAL(8), PARAMETER :: one = 1.0D0,two = 2.0d0, &
                             & zero = 0.0d0,half = 0.5d0
    REAL(8), PARAMETER :: aa = 0.5d0, bb = 0.5d0, cc = 1000.0d0, dd = 100.0d0
    REAL(8), PARAMETER :: length = 1.0D0
    REAL(8), DIMENSION(neq) :: b,phi1D
    REAL(8), DIMENSION(neq,neq) :: a
    REAL(8), DIMENSION(n) :: x,phi_bot,phi_top
    REAL(8), DIMENSION(m) :: y,phi_left,phi_right
    REAL(8), DIMENSION(n,m) :: S,phi2D,phian
    REAL(8) :: dx,dx2,dy,dy2
    REAL(4) :: ta(2),time_start,time_end
    
! Open output file
    OPEN(unit=10, file="Example3_1.out", status = "unknown")

    dx = length/(n-1); dy = length/(m-1)
    dx2 = dx*dx; dy2 = dy*dy
    
    DO i = 1,n
      x(i) = (i-1)*dx
      phi_bot(i) = cc*(bb*bb*sinh(-bb)+(x(i)-aa)*(x(i)-aa)*sinh(x(i)-aa))+dd
      phi_top(i) = cc*((one-bb)*(one-bb)*sinh(one-bb)+(x(i)-aa)*(x(i)-aa)*sinh(x(i)-aa))+dd
    ENDDO
    DO j = 1,m
      y(j) = (j-1)*dy
      phi_left(j) = cc*(aa*aa*sinh(-aa)+(y(j)-bb)*(y(j)-bb)*sinh(y(j)-bb))+dd
      phi_right(j) = cc*((one-aa)*(one-aa)*sinh(one-aa)+(y(j)-bb)*(y(j)-bb)*sinh(y(j)-bb))+dd
    ENDDO
    DO i = 1,n
      DO j = 1,m
        S(i,j) = cc*(2.0d0*sinh(x(i)-aa)+4.0d0*(x(i)-aa)*cosh(x(i)-aa)+(x(i)-aa)*(x(i)-aa)*sinh(x(i)-aa) + &
                 2.0d0*sinh(y(j)-bb)+4.0d0*(y(j)-bb)*cosh(y(j)-bb)+(y(j)-bb)*(y(j)-bb)*sinh(y(j)-bb))
      ENDDO
    ENDDO


    a(:,:) = zero; b(:) = zero

! Set up coefficient matrix and source vector

! Boundary conditions
    i = 1
    DO j = 2,m
      k = (j-1)*n+i
      a(k,k) = one
      b(k) = phi_left(j)
    ENDDO
    i = n
    DO j = 2,m
      k = (j-1)*n+i
      a(k,k) = one
      b(k) = phi_right(j)
    ENDDO
    j = 1
    DO i = 1,n
      k = (j-1)*n+i
      a(k,k) = one
      b(k) = phi_bot(i)
    ENDDO
    j = m
    DO i = 2,n-1
      k = (j-1)*n+i
      a(k,k) = one
      b(k) = phi_top(i)
    ENDDO
    DO i = 2,n-1
      DO j = 2,m-1
        k = (j-1)*n+i
        b(k) = -S(i,j)
        a(k,k-1) = -one/dx2
        a(k,k+1) = -one/dx2
        a(k,k-N) = -one/dy2
        a(k,k+N) = -one/dy2
        a(k,k) = two/dx2+two/dy2
      ENDDO
    ENDDO
     
! Solve system of equations using the Gaussian Elimination
    time_start = etime(ta)

    CALL gauss_elim(neq,a,b,phi1D)

    time_end = etime(ta)

! Print results
    WRITE(10,*) 'VARIABLES = "X", "Y", "PHIAN", "PHI", "Error"'
    WRITE(10,*) 'ZONE I=', n, ', J=', m, ', DATAPACKING=POINT'
    DO j = 1,m
      DO i = 1,n
         k = (j-1)*n+i
         phi2D(i,j) = phi1D(k)
         phian(i,j) = cc*(((x(i)-aa)**2)*sinh(x(i)-aa) + ((y(j)-bb)**2)*sinh(y(j)-bb))+dd
         WRITE(*,10) x(i),y(j),phian(i,j),phi2D(i,j),phian(i,j)-phi2D(i,j)
         WRITE(10,10) x(i),y(j),phian(i,j),phi2D(i,j),phian(i,j)-phi2D(i,j)
      ENDDO
    ENDDO

    print *, "time taken = ", time_end-time_start

 10 format(5(1x,F13.6))

    CLOSE(unit=10)

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

!-----------------------------------------------------------------------
  SUBROUTINE gauss_elim(n,a,b,x)
!-----------------------------------------------------------------------
!
! Purpose: Solves a set of linear algebraic equations.
!          Using Gaussian Elimination
!
!          n = number of equations (input)
!          a = nxn coefficient matrix (input)
!          b = right hand vector (input)
!          x = solution vector (output
!
!-----------------------------------------------------------------------

   IMPLICIT NONE
   INTEGER(4) :: i,j,k
   INTEGER(4) , INTENT(IN) :: n
   REAL(8) :: sum1,xmult
   REAL(8) :: b(n),a(n,n)
   REAL(8), INTENT(OUT) :: x(n)
!-----------------------------------------------------------------------

   DO k = 1,n-1
     DO i = k+1,n
       xmult = a(i,k)/a(k,k)
       DO j=k+1,n
         a(i,j) = a(i,j)-xmult*a(k,j)
       ENDDO
       a(i,k) = xmult
       b(i) = b(i) - xmult*b(k)
     ENDDO
   ENDDO
   x(n) = b(n)/a(n,n)
   DO i = n-1,1,-1
      sum1 = b(i)
      DO j = i+1,n
         sum1 = sum1 - a(i,j)*x(j)
      ENDDO
      x(i) = sum1/a(i,i)
   ENDDO

   END SUBROUTINE gauss_elim

