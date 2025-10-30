! Start of Main Program
! This is a program for finding the power polynomial
! fit to the function f(x) = 1/x.
! Exercise 8.2 (a)

    program main

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n = 10, m = 1001
    INTEGER(4) :: i,j,iter
    REAL(8), PARAMETER :: one = 1.0D0,two = 2.0d0, &
                             & zero = 0.0d0,half = 0.5d0
    REAL(8), DIMENSION(n) :: x,y,a
    REAL(8), DIMENSION(n,n) :: coef
    REAL(8) :: xx,dx,yexact,ypoly,sump
    
! Open output file
    OPEN(unit=10, file="Exercise8_2_poly.out", status = "unknown")

! Input data
    DO i = 1,n
      x(i) = one*i
      y(i) = one/x(i)
    ENDDO

! Set up matrix to find coeffs
    DO i = 1,n
      DO j = 1,n
        coef(i,j) = x(i)**(j-1)
      ENDDO
    ENDDO

    CALL gauss_elim(n,coef,y,a)
       
! Print results
    dx = (10.5d0-0.5d0)/(m-1)
    DO j = 1,m
      xx = 0.5d0 + (j-1)*dx
      yexact = one/xx
      sump = zero
      DO i = 1,n
        sump = sump + a(i)*(xx**(i-1))
      ENDDO
      ypoly = sump
      WRITE(10,10) xx,yexact,ypoly,100.0d0*(yexact-ypoly)/yexact
      WRITE(*,10) xx,yexact,ypoly,100.0d0*(yexact-ypoly)/yexact
    ENDDO

 10 format(F13.4,1x,2(F14.8,1x),F14.6)

    CLOSE(unit=10)

    END Program main

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
