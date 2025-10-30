! Start of Main Program
! This is a program to solve the one-dimensional advection-diffusion
! equation subject to two Dirichlet Boundary conditions using the
! finite-volume method.
! First Order Upwind Scheme
! Example 6.2

    program main

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: nc = 20
    INTEGER(4) :: iter,n_eq,i
    REAL(8), PARAMETER :: one=1.0D0,two=2.0d0,pi = 3.1415927d0, &
                             & zero= 0.0d0,phi_left=1.0D0,phi_right=0.0D0, &
                             & length=1.0D0,st=1.01d0
    REAL(8), DIMENSION(nc) :: a,b,c,d,s,x,phi
    REAL(8) :: phi_an,dx,acon,bcon,Pe
    
! Open output file
    OPEN(unit=10, file="Example6_2_1stUW.out", status = "unknown")

    Pe = 1.0d4
    
    dx = length/nc
    x(1) = dx/2.0d0
    DO i = 2,nc
      x(i) = x(i-1)+dx
    ENDDO

! Set up coefficient matrix and source vector
    d(:)=zero; a(:)=zero; b(:)=zero; c(:)=zero
    d(1) = 12.0d0+3.0d0*Pe*dx
    c(1) = -4.0d0
    b(1) = (8.0d0+3.0d0*Pe*dx)*phi_left
    DO i = 2,nc-1
      b(i) = zero
      a(i-1) = -(1.0d0+Pe*dx)
      c(i) = -1.0d0
      d(i) = 2.0d0+Pe*dx
    ENDDO
    d(nc) = 12.0d0
    a(nc-1) = -(4.0d0+3.0d0*Pe*dx)
    b(nc) = (8.0d0 - 3.0d0*Pe*dx)*phi_right

! Solve system of equations using the Thomas Algorithm
    CALL TRI(nc,a,d,c,b,phi)

! Print results
    WRITE(10,11) zero,phi_left,phi_left,zero
    WRITE(*,11) zero,phi_left,phi_left,zero
    DO i = 1,nc
      acon = (phi_left*exp(Pe)-phi_right)/(exp(Pe)-one)
      bcon = -(phi_left-phi_right)/(exp(Pe)-one)
      phi_an = acon + bcon*exp(Pe*x(i))
      WRITE(10,11) x(i),phi(i), phi_an,phi_an-phi(i)
      WRITE(*,11) x(i),phi(i), phi_an,phi_an-phi(i)
    ENDDO
    WRITE(10,11) one,phi_right,phi_right,zero
    WRITE(*,11) one,phi_right,phi_right,zero

 11 FORMAT(4(1x,E14.6))

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
