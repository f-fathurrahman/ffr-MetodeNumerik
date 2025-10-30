! Start of Main Program
! This is a program to solve a 1D diffusion equation subject
! to two Dirichlet Boundary conditions using finite difference
! Exercise 2.4

    program main

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n_nodes = 6
    INTEGER(4) :: iter,n_eq,i
    REAL(8), PARAMETER :: one=1.0D0,two=2.0d0,pi = 3.1415927d0, &
                             & zero= 0.0d0,phi_left=0.0D0,phi_right=1.0D0, &
                             & length=1.0D0
    REAL(8), DIMENSION(n_nodes) :: a,b,c,d,x,phi,phi_an
    REAL(8) :: dx,dx2,err_i
    
! Open output file
    OPEN(unit=10, file="Exercise_2_4.out", status = "unknown")

    dx = length/(n_nodes-1)
    dx2 = dx*dx
    DO i = 1,n_nodes
      x(i) = (i-1)*dx
    ENDDO

! Set up coefficient matrix and source vector
    DO i = 2,n_nodes-1
      b(i) = exp(-x(i))
      a(i-1) = one/dx2
      c(i) = one/dx2
      d(i) = -two/dx2
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

! Print results
    DO i = 2,n_nodes-1
      phi_an(i) = exp(-x(i)) + (two-exp(-one))*x(i) - one 
      err_i = 100.0d0*(phi_an(i)-phi(i))/phi_an(i)
      WRITE(10,10) x(i),phi(i), phi_an(i),err_i
    ENDDO

 10 format(4(1x,F14.6))

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
