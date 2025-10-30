! Start of Main Program
! Newton's Method
! Example 8.9

    program main

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n = 2, max_iter = 100
    INTEGER(4) :: i,iter
    REAL(8), PARAMETER :: one = 1.0D0,two = 2.0d0, &
                             & zero = 0.0d0,half = 0.5d0, &
                             & tol = 1.0D-6, omega = 0.5d0
    REAL(8), DIMENSION(n,n) :: jac
    REAL(8), DIMENSION(n) :: f,phi,dphi
    REAL(8) :: sumr,res
    
! Open output file
    OPEN(unit=10, file="Example8_9.out", status = "unknown")

    phi(1) = half
    phi(2) = half

    DO iter = 1,max_iter

      f(1) = phi(1)*phi(1)*phi(2)-exp((one-phi(1))*phi(2))-one
      f(2) = phi(1)*phi(1)*phi(1)*phi(2)+phi(2)*phi(2)*phi(1)-6.0d0
      jac(1,1) = two*phi(1)*phi(2)+phi(2)*exp((one-phi(1))*phi(2))
      jac(1,2) = phi(1)*phi(1)-(one-phi(1))*exp((one-phi(1))*phi(2))
      jac(2,1) = 3.0d0*phi(1)*phi(1)*phi(2)+phi(2)*phi(2)
      jac(2,2) = phi(1)*phi(1)*phi(1)+two*phi(1)*phi(2)

      f(:) = -f(:)

      CALL gauss_elim(n,jac,f,dphi)

      phi(:) = phi(:) + omega*dphi(:)

      sumr = zero
      DO i = 1,n
        sumr = sumr + f(i)*f(i)
      ENDDO
      res = SQRT(MAX(zero,sumr))

      WRITE(*,10) iter, res, phi(1), phi(2)
      WRITE(10,10) iter, res, phi(1), phi(2)

      IF(res < tol) EXIT

    ENDDO

    CLOSE(unit=10)

 10 format(I6,3(F14.7))

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
