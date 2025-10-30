! Start of Main Program
! This is a program to solve a 1D diffusion equation subject
! to two Dirichlet Boundary conditions using finite difference
! This problem has a non-linear source.
! Exercise 3.7

! method = 1 => source term is not linearized
! method = 2 => source term is linearized

    program main

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n_nodes = 41 , max_iter = 100, method = 2
    INTEGER(4) :: iter,n_eq,i
    REAL(8), PARAMETER :: one=1.0D0,two=2.0d0,pi = 3.1415927d0, &
                             & zero= 0.0d0,phi_left=0.0D0,phi_right=1.0D0, &
                             & length=1.0D0, tiny=1.0D-40, tol=1.0D-10, &
                             & con = 4.0d0
    REAL(8), DIMENSION(n_nodes) :: a,b,c,d,x,phi,phi_an
    REAL(8) :: dx,dx2,err_i,sumres,res,fde

    OPEN(unit=10, file="Exercise3_7.out", status = "unknown")
    OPEN(unit=20, file="Exercise3_7.rsl", status = "unknown")

    dx = length/(n_nodes-1)
    dx2 = dx*dx
    DO i = 1,n_nodes
      x(i) = (i-1)*dx
    ENDDO

! Start of iteration loop
    phi(:) = zero   ! Initial Guess for phi

    DO iter = 1,max_iter

! Set up coefficient matrix and source vector
      DO i = 2,n_nodes-1
        b(i) = exp(con*phi(i))
        IF(method == 2) b(i) = b(i)-con*phi(i)*exp(con*phi(i))
        a(i-1) = one/dx2
        c(i) = one/dx2
        d(i) = -two/dx2
        IF(method == 2) d(i) = d(i)-con*exp(con*phi(i))
      ENDDO

! Boundary conditions
      d(1) = one
      c(1) = zero
      b(1) = phi_left
      d(n_nodes) = one
      a(n_nodes-1) = zero
      c(n_nodes) = zero
      b(n_nodes) = phi_right

! Solve system of equations using the Thomas Algorithm
      CALL TRI(n_nodes,a,d,c,b,phi)

!     Calculate L2NORM
      sumres = zero
      DO i = 2,n_nodes-1
        fde = (phi(i+1)-two*phi(i)+phi(i-1))/dx2-exp(con*phi(i))
        sumres = sumres + fde*fde
      ENDDO
      res = sqrt(MAX(tiny,sumres))

      WRITE(*,*) iter, res
      WRITE(20,*) iter, res
      IF(res < tol) EXIT    ! Convergence check

    ENDDO ! iteration loop

! Print results
    DO i = 1,n_nodes
      WRITE(10,*) x(i),phi(i)
    ENDDO

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
