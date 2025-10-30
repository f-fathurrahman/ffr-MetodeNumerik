! Start of Main Program
! This is a program to solve an initial value
! problem in two independent variables using the
! Backward Euler method and Newton iterations
! Exercise 8.6

    program main

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n = 2 , max_iter = 100, max_times = 100
    INTEGER(4) :: iter,i,it
    REAL(8), PARAMETER :: one=1.0D0,two=2.0d0,zero= 0.0d0, &
                             & tol=1.0D-10, omega = 1.0d0
    REAL(8), DIMENSION(n) :: f,phi,phio,rhs,dphi
    REAL(8), DIMENSION(n,n) :: jac
    REAL(8) :: err_i,sumres,res,fde,dt,t
    
! Open output file
    OPEN(unit=10, file="Exercise8_6.out", status = "unknown")

    dt = 0.01d0
    phio(:) = one ! Initial condition
    phi(:) = phio(:)   ! Initial Guess for phi

! Start of time advancement
    t = zero

    DO it = 1,max_times

      t = t + dt

! Start of Newton Iteration

      DO iter = 1,max_iter

! Calculate function
        f(1) = (phi(1) - phio(1))/dt - 10.0d0*(phi(1)**2)*exp(-one/phi(2))
        f(2) = (phi(2) - phio(2))/dt + phi(2) - one - phi(1)*exp(-one/phi(2))
        
        rhs(:) = -f(:) ! Right hand side of equations

! Calculate Jacobian (Coefficient Matrix)
        jac(1,1) = one/dt - 20.0d0*phi(1)*exp(-one/phi(2))
        jac(1,2) = - 10.0d0*(phi(1)**2)*exp(-one/phi(2))/(phi(2)**2)
        jac(2,1) = - exp(-one/phi(2))
        jac(2,2) = one/dt + one - phi(1)*exp(-one/phi(2))/(phi(2)**2)

! Solve system of equations using Gaussian Elimination
        CALL gauss_elim(n,jac,rhs,dphi)

! Update Solution
        phi(:) = phi(:) + omega*dphi(:)

! Calculate L2NORM
        sumres = zero
        DO i = 1,n
          sumres = sumres + f(i)*f(i)
        ENDDO
        res = sqrt(MAX(zero,sumres))

        print *, iter, res
        IF(res < tol) EXIT    ! Convergence check

      ENDDO ! iteration loop

      phio(:) = phi(:) ! Reset initial condition

! Print results
      WRITE(10,*) t, (phi(i), i = 1,n)
      WRITE(*,*) t, (phi(i), i = 1,n)

    ENDDO ! Time loop


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
