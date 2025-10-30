! Start of Main Program
! Example 6.1
! Finite Volume Method

    program main

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n_cells = 8
    INTEGER(4) :: i
    REAL(8), PARAMETER :: one=1.0D0,two=2.0d0,zero= 0.0d0
    REAL(8), PARAMETER :: phi_left=0.0d0, phi_right=0.0D0
    REAL(8), DIMENSION(n_cells) :: a,d,c,b,s,phi,x,phian
    REAL(8) :: dx,dx2,length,J_l,J_r,c1,c2,sphi,imb
    REAL(8) :: J_an_l,J_an_r,sphi_an,imb_an
    
! Open output file
    OPEN(unit=10, file="Example6_1_fvm.out", status = "unknown")

    length = one
    dx = length/(n_cells)
    dx2 = dx*dx
    x(1) = dx/two
    s(1) = cos(x(1))
    DO i = 2,n_cells
      x(i) = x(1)+(i-1)*dx
      s(i) = cos(x(i))
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
! Phi
    c2 = one
    c1 = (-c2+cos(length))/length
    WRITE(*,10) zero,phi_left,phi_left,zero
    WRITE(10,10) zero,phi_left,phi_left,zero
    DO i = 1,n_cells
      phian(i) = -cos(x(i)) + c1*x(i) + c2
      WRITE(*,10) x(i),phian(i),phi(i),phian(i)-phi(i)
      WRITE(10,10) x(i),phian(i),phi(i),phian(i)-phi(i)
    ENDDO
    WRITE(*,10) length,phi_right,phi_right,zero
    WRITE(10,10) length,phi_right,phi_right,zero
! Fluxes
    J_l = -(9.0d0*phi(1)-phi(2)-8.0d0*phi_left)/(3.0d0*dx)
    J_r = (9.0d0*phi(n_cells)-phi(n_cells-1)-8.0d0*phi_right)/(3.0d0*dx)
    sphi = zero
    DO i = 1,n_cells
      sphi = sphi + s(i)*dx
    ENDDO
    sphi = -sphi
    imb = J_l-J_r+sphi

    WRITE(*,*) "  Left          ","Right         ","Source          ","Imbalance   "
    WRITE(10,*)"  Left          ","Right         ","Source          ","Imbalance   "
    WRITE(*,10) J_l,J_r,sphi,imb
    WRITE(10,10) J_l,J_r,sphi,imb

    J_an_l = -(c1)
    J_an_r = -(sin(length)+c1)
    sphi_an = -sin(length)
    imb_an = J_an_l-J_an_r+sphi_an
    WRITE(*,*) "Analytical Solution"
    WRITE(*,10) J_an_l, J_an_r, sphi_an,imb_an
    WRITE(10,*) "Analytical Solution"
    WRITE(10,10) J_an_l, J_an_r, sphi_an,imb_an

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


