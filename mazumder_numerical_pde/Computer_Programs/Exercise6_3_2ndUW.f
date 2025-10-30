! Start of Main Program
! This is a program to solve the one-dimensional advection-diffusion
! equation subject to two Dirichlet Boundary conditions using the
! finite-volume method.
! Second Order UDS
! Exercise 6.3 (ii)

    program main

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: nc = 160
    INTEGER(4) :: iter,n_eq,i
    REAL(8), PARAMETER :: one=1.0D0,two=2.0d0,pi = 3.1415927d0, &
                             & zero= 0.0d0,phi_left=0.0D0,phi_right=1.0D0, &
                             & length=1.0D0
    REAL(8), DIMENSION(nc) :: b,x,phi
    REAL(8), DIMENSION(nc,nc) :: a
    REAL(8) :: phi_an,dx,acon,bcon,Pe
    
! Open output file
    OPEN(unit=10, file="Exercise6_3_2ndUW.out", status = "unknown")

    Pe = 1.0d1
    
    dx = length/nc
    x(1) = dx/2.0d0
    DO i = 2,nc
      x(i) = x(i-1)+dx
    ENDDO

! Set up coefficient matrix and source vector
    a(:,:)=zero; b(:) = zero
    a(1,1) = 12.0d0+6.0d0*Pe*dx
    a(1,2) = -4.0d0
    b(1) = (8.0d0+6.0d0*Pe*dx)*phi_left
    a(2,1) = -(2.0d0+5.0d0*Pe*dx)
    a(2,2) = (4.0d0+3.0d0*Pe*dx)
    a(2,3) = -2.0d0
    b(2) = -2.0d0*Pe*dx*phi_left
    DO i = 3,nc-1
      a(i,i-1) = -(2.0d0+4.0d0*Pe*dx)
      a(i,i) = 4.0d0+3.0d0*Pe*dx
      a(i,i+1) = -2.0d0
      a(i,i-2) = Pe*dx
      b(i) = zero
    ENDDO
    a(nc,nc) = 24.0d0
    a(nc,nc-1) = -(8.0d0+9.0d0*Pe*dx)
    a(nc,nc-2) = 3.0d0*Pe*dx
    b(nc) = (16.0d0 - 6.0d0*Pe*dx)*phi_right

! Solve system of equations using the Thomas Algorithm
    CALL gauss_elim(nc,a,b,phi)

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

