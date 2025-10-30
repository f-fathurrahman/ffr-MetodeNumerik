! Start of Main Program
! This is a program to solve a 1D linear ODE using the FVM.
! Exercise 6.1

    program main

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n_cells = 20
    INTEGER(4) :: i
    REAL(8), PARAMETER :: one=1.0D0,two=2.0d0,zero= 0.0d0
    REAL(8), PARAMETER :: phi_left=0.0d0, phi_right=1.0D0
    REAL(8), DIMENSION(n_cells) :: a,d,c,b,s,phi,x,phian
    REAL(8) :: dx,dx2,length
    
! Open output file
    OPEN(unit=10, file="Exercise6_1.out", status = "unknown")

    length = one
    dx = length/(n_cells)
    dx2 = dx*dx
    x(1) = dx/two
    s(1) = exp(x(1))
    DO i = 2,n_cells
      x(i) = x(1)+(i-1)*dx
      s(i) = exp(x(i))
    ENDDO

! Set up coefficient matrix and source vector

    a(:) = zero; b(:) = zero; c(:) = zero; d(:) = zero; phi(:) = zero

! Boundary conditions

    d(1) = 4.0d0/dx2
    c(1) = -4.0d0/(3.0d0*dx2)
    b(1) = -s(1)+8.0d0*phi_left/(3.0d0*dx2)
    DO i = 2,n_cells-1
      b(i) = -s(i)
      d(i) = two/dx2
      c(i) = -one/dx2
      a(i-1) = -one/dx2
    ENDDO
    d(n_cells) = 4.0d0/dx2
    a(n_cells-1) = -4.0d0/(3.0d0*dx2)
    b(n_cells) = -s(n_cells)+8.0d0*phi_right/(3.0d0*dx2)

! Solve system of equations using the TDMA

    CALL TRI(n_cells,a,d,c,b,phi)

! Print results
    WRITE(*,10) zero,phi_left,phi_left,zero
    WRITE(10,10) zero,phi_left,phi_left,zero
    DO i = 1,n_cells
      phian(i) = exp(x(i))+(two-exp(one))*x(i)-one
      WRITE(*,10) x(i),phian(i),phi(i),phian(i)-phi(i)
      WRITE(10,10) x(i),phian(i),phi(i),phian(i)-phi(i)
    ENDDO
    WRITE(*,10) length,phi_right,phi_right,zero
    WRITE(10,10) length,phi_right,phi_right,zero

 10 format(4(F14.6))

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


